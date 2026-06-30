import tempfile
import warnings
import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as smf
from volgrids._vendors import freyacli as fy
from volgrids._vendors import molsimple as ms

# //////////////////////////////////////////////////////////////////////////////
class MoleculeManager:
    """Manages the molecular system, including the structure, (optional) trajectory, enclosing box, chemtable and selection queries."""

    # --------------------------------------------------------------------------
    def __init__(self, path_struct: Path, path_traj: Path = None, box: vg.Box = None):
        import MDAnalysis as mda

        self.molname  : str                # name of the molecule
        self.chemtable: smf.ParserChemTable
        self.particles: ms.ParticleGroup

        self.frame     : int|None           # current frame number (if trajectory is used)
        self.nframes   : int                # total number of frames in the trajectory (1 if no trajectory is used)
        self._mda_universe: mda.Universe|None  # MDAnalysis Universe object for the molecular system

        self.box       : vg.Box             # current box (either enforced, sphere-based, or computed from the structure)
        self.box_common: vg.Box|None        # box that can enclose all frames of a trajectory

        self.do_traj            : bool # whether this is a trajectory or a single structure (None if no structure is provided)
        self.do_use_sphere      : bool
        self.do_enforce_boxes   : bool
        self.do_use_common_box  : bool
        self.enforce_cmap_output: bool

        self.molname = path_struct.stem

        self.do_traj = path_traj is not None
        self.do_use_sphere = len(smf.SPHERES) > 0
        self.do_enforce_boxes = len(smf.BOXES_ENFORCED) > 0
        self.do_use_common_box = self.do_traj and not vg.CFG.box_tight_traj \
            and not self.do_use_sphere and not self.do_enforce_boxes
        self.enforce_cmap_output = self.do_traj # can get updated with other criteria too (e.g. --pack flag in app_smiffer.py)

        self.chemtable = self._init_chemtable()

        self._mda_universe = vg.Utils.create_mda_universe_quiet(path_struct, path_traj)
        self.particles = ms.System(path_struct).particles # in the case of multiple models: `System.particles` is the first model

        self.frame = 0 if self.do_traj else None
        self.nframes = self._mda_universe.trajectory.n_frames if self.do_traj else 1

        if self.do_use_sphere:
            vg.SphereInfo.assert_sphere_infos(smf.SPHERES, self.nframes)

        if self.do_enforce_boxes:
            vg.BoxInfo.assert_box_infos([box.info for box in smf.BOXES_ENFORCED], self.nframes)

        self.box_common = self._gen_box_common() if self.do_use_common_box else None
        self.box = self._get_current_box() if box is None else box


    # --------------------------------------------------------------------------
    @classmethod
    def from_pqr_data(cls, pqr_data: str, box: vg.Box = None, chains: list[str] = None):
        """Adds back the `chains` information that is empty in the PQR file. `chains` should be of size (nresidues,)."""
        if not pqr_data:
            raise ValueError("Empty PQR content, aborting MoleculeManager instantiation.")

        with tempfile.NamedTemporaryFile(mode = "w+", suffix = ".pqr", delete = True) as tmp_pqr:
            tmp_pqr.write(pqr_data)
            tmp_pqr.flush()
            obj = cls(Path(tmp_pqr.name), path_traj = None, box = box)

        ### add back the chain information
        if chains is None:
            chains = ['A'] * len(obj._mda_universe.residues) # [TODO] remove MDA

        chains_per_atom = [
            chains[i]
            for i,residue in enumerate(obj._mda_universe.residues) # [TODO] remove MDA
            for _ in residue.atoms
        ]

        obj._mda_universe.add_TopologyAttr("chainIDs", chains_per_atom) # [TODO] remove MDA
        return obj


    # --------------------------------------------------------------------------
    @staticmethod
    def copy_attrs_except_universe(src: "MoleculeManager", dst: "MoleculeManager"):
        ### [TODO]: rework this method? it's easy to forget to add new attributes here when adding them to the class
        dst.molname    = src.molname
        dst.do_traj    = src.do_traj
        dst.frame      = src.frame
        dst.box        = src.box
        dst.box_common = src.box_common
        dst.chemtable  = src.chemtable
        dst.do_use_sphere       = src.do_use_sphere
        dst.do_enforce_boxes    = src.do_enforce_boxes
        dst.do_use_common_box   = src.do_use_common_box
        dst.enforce_cmap_output = src.enforce_cmap_output


    # --------------------------------------------------------------------------
    def switch_frame(self, frame_idx: int):
        self._mda_universe.trajectory[frame_idx]
        self.frame = frame_idx
        self.box = self._get_current_box()


    # --------------------------------------------------------------------------
    def get_min_coords(self): return self.box.min_coords
    def get_max_coords(self): return self.box.max_coords
    def get_resolution(self): return self.box.resolution
    def get_deltas(self):     return self.box.deltas
    def get_cog(self):        return self.box.cog
    def get_radius(self):     return self.box.radius


    # --------------------------------------------------------------------------
    def get_all_atoms(self):
        """Returns all atoms in the molecular system, without any filtering (e.g. selection query, sphere, etc.)."""
        return self._mda_universe.atoms # [TODO] remove MDA


    # --------------------------------------------------------------------------
    def get_all_queried_atoms(self, use_custom = True):
        """Returns all atoms that match the selection query, without any additional filtering (e.g. sphere)."""
        query = self.chemtable.get_selection_query(use_custom)
        atoms = self._mda_universe.select_atoms(query) # [TODO] remove MDA
        if len(atoms) == 0: warnings.warn(
            f"\n\n... The selection query '{fy.Color.blue(query)}' {fy.Color.red('did not return any atoms')}."
        )
        return atoms


    # --------------------------------------------------------------------------
    def get_relevant_queried_atoms(self, use_custom = True, extra_dist: float = 0.0):
        """Returns all atoms that match the selection query, with additional filtering (e.g. sphere)."""
        query = self.chemtable.get_selection_query(use_custom)
        if self.do_use_sphere:
            sphere = smf.SPHERES[self.frame or 0]
            query += f"and point {sphere.get_str_query(extra_dist)}"

        atoms = self._mda_universe.select_atoms(query) # [TODO] remove MDA
        if len(atoms) == 0: warnings.warn(
            f"\n\n... The selection query '{fy.Color.blue(query)}' {fy.Color.red('did not return any atoms')}."
        )
        return atoms


    # --------------------------------------------------------------------------
    def get_hydrogens(self):
        return self._mda_universe.select_atoms("not water and name H*") # [TODO] remove MDA


    # --------------------------------------------------------------------------
    def get_residue_chains(self) -> list[str]:
        """Returns a list of chain IDs for each residue in the molecular system."""
        return [
            arr[0] for arr in self._mda_universe.residues.chainIDs # [TODO] remove MDA
        ] # size: (nresidues,)


    # --------------------------------------------------------------------------
    def _get_current_box(self) -> vg.Box:
        ### Priority 1: Box to be used is specifically requested by the user (via CLI)
        if self.do_enforce_boxes:
            return smf.BOXES_ENFORCED[self.frame or 0]

        ### Priority 2: Box to be used is based on a user-provided sphere (via CLI)
        if self.do_use_sphere:
            sphere = smf.SPHERES[self.frame or 0]
            box = vg.Box(None, None, None, do_init = False)
            box.cog = np.array(sphere.get_pos())
            box.min_coords = box.cog - sphere.radius
            box.max_coords = box.cog + sphere.radius
            box.radius = sphere.radius
            box.infer_deltas_resolution()
            if vg.CFG.box_force_equilateral: box.enforce_equilateral()
            return box

        ### Priority 3: dealing with a trajectory and the user requested the box
        ### to be inferred from the structure at every frame (via config `BOX_TIGHT_TRAJ=True`)
        if not self.do_use_common_box:
            ### [TODO] MDA has to be kept here for dealing with traj --> update coords of self.particles with the current frame's coordinates
            min_coords = self._mda_universe.coord.positions.min(axis = 0)
            max_coords = self._mda_universe.coord.positions.max(axis = 0)
            return self._padded_box(min_coords, max_coords)

        ### Priority 4 (default): a common box that can enclose the structure in all frames is used (via config `BOX_TIGHT_TRAJ=False`)
        return self.box_common


    # --------------------------------------------------------------------------
    def _gen_box_common(self) -> vg.Box:
        """
        In the case of trajectories and when `vg.CFG.box_tight_traj=False`, a box that can enclose all frames is computed.
        This behavior is overridden if the user explicitly requests a box (via CLI) or a sphere (via CLI) to be used.
        """
        min_coords = np.full(3,  np.inf)
        max_coords = np.full(3, -np.inf)
        for _ in self._mda_universe.trajectory:
            positions = self._mda_universe.coord.positions
            np.minimum(min_coords, positions.min(axis = 0), out = min_coords)
            np.maximum(max_coords, positions.max(axis = 0), out = max_coords)
        self._mda_universe.trajectory[0] # rewind to frame 0
        return self._padded_box(min_coords, max_coords)


    # --------------------------------------------------------------------------
    @staticmethod
    def _padded_box(min_coords: np.ndarray, max_coords: np.ndarray) -> vg.Box:
        box = vg.Box.from_min_max(
            min_coords = min_coords - vg.CFG.box_padding,
            max_coords = max_coords + vg.CFG.box_padding,
        )
        if vg.CFG.box_force_equilateral: box.enforce_equilateral()
        return box


    # --------------------------------------------------------------------------
    def _init_chemtable(self) -> "smf.ParserChemTable":
        if smf.PATH_CHEM_LIGAND:
            return smf.ParserChemTable(smf.PATH_CHEM_LIGAND)

        folder_default_tables = vg.Utils.resolve_path_package("_data/smiffer_tables")
        chem = smf.ParserChemTable(folder_default_tables / "default.chem")

        if vg.CFG.smif_hb_only_nbase:
            ini = vg.ParserIni.from_file(folder_default_tables / "nucl_simple_hb.chem")
            chem.parse_names_hbacceptors(ini)
            chem.parse_names_hbdonors(ini)

        return chem


# //////////////////////////////////////////////////////////////////////////////
