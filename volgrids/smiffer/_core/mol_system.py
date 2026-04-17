import tempfile
import warnings
import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class MolSystemSmiffer(vg.MolSystem):
    def __init__(self, path_struct: Path, path_traj: Path = None):
        self.do_ps = sm.SPHERE is not None
        self.chemtable = sm.ParserChemTable(self._get_path_table())
        super().__init__(path_struct, path_traj)


    # --------------------------------------------------------------------------
    @classmethod
    def from_pqr_data(cls, pqr_data: str):
        with tempfile.NamedTemporaryFile(mode = "w+", suffix = ".pqr", delete = True) as tmp_pqr:
            tmp_pqr.write(pqr_data)
            tmp_pqr.flush()
            return cls(Path(tmp_pqr.name), None)


    # --------------------------------------------------------------------------
    @staticmethod
    def copy_attributes_except_system(src: "MolSystemSmiffer", dst: "MolSystemSmiffer"):
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
            f"\n\n... The selection query '{query}' did not return any atoms. "+\
            "Are you using the adequate 'prot'/'rna'/'ligand' mode?\n"
        )
        return atoms


    # --------------------------------------------------------------------------
    def _get_init_box(self) -> vg.Box: # override base class to add pocket_sphere behavior
        if self.do_ps:
            box = vg.Box(None, None, None, do_init = False)
            box.cog = np.array([sm.SPHERE.x, sm.SPHERE.y, sm.SPHERE.z])
            box.min_coords = box.cog - sm.SPHERE.radius
            box.max_coords = box.cog + sm.SPHERE.radius
            box.radius = sm.SPHERE.radius
            box.infer_deltas_resolution()
            return box

        return super()._get_init_box()


    # --------------------------------------------------------------------------
    def _get_path_table(self) -> Path:
        if sm.PATH_TABLE: return sm.PATH_TABLE

        folder_default_tables = Path("_tables")

        if sm.CURRENT_MOLTYPE == sm.MolType.PROT:
            return vg.resolve_path_package(folder_default_tables / "prot.chem")

        if sm.CURRENT_MOLTYPE == sm.MolType.RNA:
            name = "rna_simple_hb" if sm.HBONDS_ONLY_NUCLEOBASE else "rna"
            return vg.resolve_path_package(folder_default_tables / f"{name}.chem")

        raise ValueError(f"No default table for the specified molecular type '{sm.CURRENT_MOLTYPE}'. Please provide a path to a custom table.")


# //////////////////////////////////////////////////////////////////////////////
