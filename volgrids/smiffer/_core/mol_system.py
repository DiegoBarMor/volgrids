import tempfile
import warnings
import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class MolSystem:
    def __init__(self, path_struct: Path, path_traj: Path = None, box: vg.Box = None):
        import MDAnalysis as mda

        self.molname  : str                 # name of the molecule
        self.do_traj  : None | bool         # whether this is a trajectory or a single structure (None if no structure is provided)
        self.system   : None | mda.Universe # MDAnalysis Universe object for the molecular system
        self.frame    : None | int          # current frame number (if trajectory is used)
        self.box      : vg.Box
        self.do_ps    : bool
        self.chemtable: sm.ParserChemTable

        self.molname = path_struct.stem
        self.do_traj = path_traj is not None
        self.do_ps = len(sm.SPHERES) > 0
        self.do_box_enforced = sm.BOX_ENFORCED is not None
        self.chemtable = self._init_chemtable()

        self.system = vg.Utils.create_mda_universe_quiet(path_struct, path_traj)
        self.frame = 0 if self.do_traj else None

        self.box = self._get_init_box() if box is None else box

        if self.do_ps:
            nframes = self.system.trajectory.n_frames if self.do_traj else 1
            sm.SphereInfo.assert_sphere_list(sm.SPHERES, nframes)

        self.enforce_cmap_output = self.do_traj # can get updated with other criteria too (e.g. --pack flag in app_smiffer.py)


    # --------------------------------------------------------------------------
    @classmethod
    def from_pqr_data(cls, pqr_data: str, box: vg.Box = None, chains: list[str] = None):
        """Adds back the `chains` information that is empty in the PQR file. `chains` should be of size (nresidues,)."""
        if not pqr_data:
            raise ValueError("Empty PQR content, aborting MolSystem instantiation.")

        with tempfile.NamedTemporaryFile(mode = "w+", suffix = ".pqr", delete = True) as tmp_pqr:
            tmp_pqr.write(pqr_data)
            tmp_pqr.flush()
            obj = cls(Path(tmp_pqr.name), path_traj = None, box = box)

        ### add back the chain information
        if chains is None:
            chains = ['A'] * len(obj.system.residues)

        chains_per_atom = [
            chains[i]
            for i,residue in enumerate(obj.system.residues)
            for _ in residue.atoms
        ]

        obj.system.add_TopologyAttr("chainIDs", chains_per_atom)
        return obj



    # --------------------------------------------------------------------------
    @staticmethod
    def copy_attributes_except_system(src: "MolSystem", dst: "MolSystem"):
        dst.molname = src.molname
        dst.do_traj = src.do_traj
        dst.frame = src.frame
        dst.box = src.box
        dst.do_ps = src.do_ps
        dst.chemtable = src.chemtable
        dst.enforce_cmap_output = src.enforce_cmap_output


    # --------------------------------------------------------------------------
    def switch_frame(self, frame_idx: int):
        self.system.trajectory[frame_idx]
        self.frame = frame_idx


    # --------------------------------------------------------------------------
    def get_min_coords(self): return self.box.min_coords
    def get_max_coords(self): return self.box.max_coords
    def get_resolution(self): return self.box.resolution
    def get_deltas(self):     return self.box.deltas
    def get_cog(self):        return self.box.cog
    def get_radius(self):     return self.box.radius


    # --------------------------------------------------------------------------
    def get_all_atoms(self, use_custom = True):
        query = self.chemtable.get_selection_query(use_custom)
        atoms = self.system.select_atoms(query)
        if len(atoms) == 0: warnings.warn(
            f"\n\n... The selection query '{fy.Color.blue(query)}' {fy.Color.red('did not return any atoms')}."
        )
        return atoms


    # --------------------------------------------------------------------------
    def get_relevant_atoms(self, use_custom = True, extra_dist: float = 0.0):
        query = self.chemtable.get_selection_query(use_custom)
        if self.do_ps:
            sphere = self._get_current_sphere()
            query += f"and point {sphere.get_str_query(extra_dist)}"

        atoms = self.system.select_atoms(query)
        if len(atoms) == 0: warnings.warn(
            f"\n\n... The selection query '{fy.Color.blue(query)}' {fy.Color.red('did not return any atoms')}."
        )
        return atoms


    # --------------------------------------------------------------------------
    def _get_init_box(self) -> vg.Box:
        if self.do_box_enforced: return sm.BOX_ENFORCED

        if self.do_ps:
            sphere = self._get_current_sphere()
            box = vg.Box(None, None, None, do_init = False)
            box.cog = np.array(sphere.get_pos())
            box.min_coords = box.cog - sphere.radius
            box.max_coords = box.cog + sphere.radius
            box.radius = sphere.radius
            box.infer_deltas_resolution()
            if vg.ENSURE_EQUILATERAL: box.enforce_equilateral()
            return box

        if self.do_traj:
            min_coords = np.full(3,  np.inf)
            max_coords = np.full(3, -np.inf)
            for _ in self.system.trajectory:
                positions = self.system.coord.positions
                np.minimum(min_coords, positions.min(axis = 0), out = min_coords)
                np.maximum(max_coords, positions.max(axis = 0), out = max_coords)
            self.system.trajectory[0] # rewind to frame 0
        else:
            min_coords = np.min(self.system.coord.positions, axis = 0)
            max_coords = np.max(self.system.coord.positions, axis = 0)

        box = vg.Box.from_min_max(
            min_coords = min_coords - vg.EXTRA_BOX_SIZE,
            max_coords = max_coords + vg.EXTRA_BOX_SIZE,
        )
        if vg.ENSURE_EQUILATERAL: box.enforce_equilateral()
        return box


    # --------------------------------------------------------------------------
    def _init_chemtable(self) -> "sm.ParserChemTable":
        if sm.PATH_CHEM_LIGAND:
            return sm.ParserChemTable(sm.PATH_CHEM_LIGAND)

        folder_default_tables = vg.Utils.resolve_path_package("_data/smiffer_tables")
        chem = sm.ParserChemTable(folder_default_tables / "default.chem")

        if sm.HBONDS_ONLY_NUCLEOBASE:
            ini = vg.ParserIni.from_file(folder_default_tables / "nucl_simple_hb.chem")
            chem.parse_names_hbacceptors(ini)
            chem.parse_names_hbdonors(ini)

        return chem


    # --------------------------------------------------------------------------
    def _get_current_sphere(self) -> "sm.SphereInfo":
        if not self.do_ps: raise ValueError("No pocket sphere info available.")
        if self.frame is None: return sm.SPHERES[0]
        return sm.SPHERES[self.frame]



# //////////////////////////////////////////////////////////////////////////////
