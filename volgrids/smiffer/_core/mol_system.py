import tempfile
import warnings
import numpy as np
import MDAnalysis as mda
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class MolSystem:
    def __init__(self, path_struct: Path, path_traj: Path = None, box: vg.Box = None):
        self.molname  : str                 # name of the molecule
        self.do_traj  : None | bool         # whether this is a trajectory or a single structure (None if no structure is provided)
        self.system   : None | mda.Universe # MDAnalysis Universe object for the molecular system
        self.frame    : None | int          # current frame number (if trajectory is used)
        self.box      : vg.Box
        self.do_ps    : bool
        self.chemtable: sm.ParserChemTable

        self.molname = path_struct.stem
        self.do_traj = path_traj is not None
        self.do_ps = sm.SPHERE is not None
        self.chemtable = self._init_chemtable()

        self.system = vg.Utils.create_mda_universe_quiet(path_struct, path_traj)
        self.frame = 0 if self.do_traj else None

        self.box = self._get_init_box() if box is None else box


    # --------------------------------------------------------------------------
    @classmethod
    def from_pqr_data(cls, pqr_data: str, box: vg.Box = None):
        if not pqr_data:
            raise ValueError("Empty PQR content, aborting MolSystem instantiation.")

        with tempfile.NamedTemporaryFile(mode = "w+", suffix = ".pqr", delete = True) as tmp_pqr:
            tmp_pqr.write(pqr_data)
            tmp_pqr.flush()
            return cls(Path(tmp_pqr.name), path_traj = None, box = box)


    # --------------------------------------------------------------------------
    @staticmethod
    def copy_attributes_except_system(src: "MolSystem", dst: "MolSystem"):
        dst.molname = src.molname
        dst.do_traj = src.do_traj
        dst.frame = src.frame
        dst.box = src.box
        dst.do_ps = src.do_ps
        dst.chemtable = src.chemtable


    # --------------------------------------------------------------------------
    def get_min_coords(self): return self.box.min_coords
    def get_max_coords(self): return self.box.max_coords
    def get_resolution(self): return self.box.resolution
    def get_deltas(self):     return self.box.deltas
    def get_cog(self):        return self.box.cog
    def get_radius(self):     return self.box.radius


    # --------------------------------------------------------------------------
    def get_relevant_atoms(self, use_custom = True, extra_dist: float = 0.0):
        query = self.chemtable.get_selection_query(use_custom)
        if self.do_ps: query += f"and point {sm.SPHERE.get_str_query(extra_dist)}"

        atoms = self.system.select_atoms(query)
        if len(atoms) == 0: warnings.warn(
            f"\n\n... The selection query '{fy.Color.blue(query)}' {fy.Color.red('did not return any atoms')}."
        )
        return atoms


    # --------------------------------------------------------------------------
    def _get_init_box(self) -> vg.Box:
        if self.do_ps:
            box = vg.Box(None, None, None, do_init = False)
            box.cog = np.array([sm.SPHERE.x, sm.SPHERE.y, sm.SPHERE.z])
            box.min_coords = box.cog - sm.SPHERE.radius
            box.max_coords = box.cog + sm.SPHERE.radius
            box.radius = sm.SPHERE.radius
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
            ini = vg.ParserIni.from_file(folder_default_tables / "rna_simple_hb.chem")
            chem.parse_names_hbacceptors(ini)
            chem.parse_names_hbdonors(ini)

        return chem


# //////////////////////////////////////////////////////////////////////////////
