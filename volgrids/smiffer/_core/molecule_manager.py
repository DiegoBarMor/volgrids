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

        self.atoms_all: ms.ParticleGroup # all atoms in the molecular system, without any filtering (e.g. selection query, sphere, etc.)
        self.atoms_filter_trim: ms.ParticleGroup # atoms from the requested resname (from the chemtable) and with no hydrogens, to be used by Trimmer
        self.atoms_filter_smif: ms.ParticleGroup # same as `atoms_filter_trim` but further filtered to only include atoms part of the optional explicit filtering of custom residues (from the CLI)

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
        self.init_atoms(path_struct)

        self.frame = 0 if self.do_traj else None
        self.nframes = self._mda_universe.trajectory.n_frames if self.do_traj else 1

        if self.do_use_sphere:
            vg.SphereInfo.assert_sphere_infos(smf.SPHERES, self.nframes)

        if self.do_enforce_boxes:
            vg.BoxInfo.assert_box_infos([box.info for box in smf.BOXES_ENFORCED], self.nframes)

        self.box_common = self._gen_box_common() if self.do_use_common_box else None
        self.box = self._get_current_box() if box is None else box


    # --------------------------------------------------------------------------
    def init_atoms(self, path_struct: Path, chains: list[str] = None):
        """Can add back the `chains` information that is empty in PQR files. `chains` should be of size (nresidues,)."""
        self.atoms_all = ms.System.read_pdb(path_struct).particles # in the case of multiple models: `System.atoms_all` is the first model
        smf.ResnameStandard.standardize_particle_group(self.atoms_all)
        self.atoms_filter_trim = self._init_filter_trim()

        if chains is not None:
            residues = self.atoms_all.split_residues()
            chains_per_atom = [
                chains[i] for i,residue in enumerate(residues) for _ in residue
            ]
            self.atoms_all.set_chainids(chains_per_atom)

        self.atoms_filter_smif = self._init_filter_smif()


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
    def get_atoms_insphere(self, extra_dist: float = 0.0) -> "ms.ParticleGroup":
        """
        Same as `atoms_filter_smif` but further filtered to only include atoms inside an (optional) sphere.
        It must be calculated at every frame if the sphere is moving (e.g. when using a trajectory).
        """
        atoms = self.atoms_filter_smif
        if self.do_use_sphere:
            sphere = smf.SPHERES[self.frame or 0]
            atoms = sphere.filter_particles(atoms, extra_dist)
        return atoms


    # --------------------------------------------------------------------------
    def get_hydrogens(self):
        return self._mda_universe.select_atoms("not water and name H*") # [TODO] remove MDA


    # --------------------------------------------------------------------------
    def get_residue_chains(self) -> list[str]:
        """Returns a list of chain IDs for each residue in the molecular system, size: (nresidues,)"""
        return [residue[0].chainid for residue in self.atoms_all.split_residues()]


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
            ### [TODO] MDA has to be kept here for dealing with traj --> update coords of self.atoms_all with the current frame's coordinates
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


    # --------------------------------------------------------------------------
    def _init_filter_trim(self) -> "ms.ParticleGroup":
        atoms = self.atoms_all\
            .select_resname(*self.chemtable.resnames)\
            .select_non_hydrogens()

        if len(atoms) == 0: warnings.warn(
            f"\n\n... {fy.Color.red('Empty atom selection')}."
        )
        return atoms


    # --------------------------------------------------------------------------
    def _init_filter_smif(self):
        atoms = self.atoms_filter_trim

        if smf.CUSTOM_RESIDUES:
            chainresids = (ms.ChainResid.from_dotstr(s) for s in smf.CUSTOM_RESIDUES.split())
            atoms = atoms.select_chain_resid(*chainresids)

        if len(atoms) == 0: warnings.warn(
            f"\n\n... {fy.Color.red('Empty atom selection')}."
        )
        return atoms


# //////////////////////////////////////////////////////////////////////////////
