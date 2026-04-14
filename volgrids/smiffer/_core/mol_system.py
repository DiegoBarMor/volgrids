import tempfile
import numpy as np
from pathlib import Path
from enum import Enum, auto

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class MolType(Enum):
    NONE = auto()
    PROT = auto()
    RNA = auto()
    LIGAND = auto()

    # --------------------------------------------------------------------------
    def is_none(self):   return self == MolType.NONE
    def is_prot(self):   return self == MolType.PROT
    def is_rna(self):    return self == MolType.RNA
    def is_ligand(self): return self == MolType.LIGAND


# //////////////////////////////////////////////////////////////////////////////
class MolSystemSmiffer(vg.MolSystem):
    def __init__(self, path_struct: Path, path_traj: Path = None):
        self.do_ps = sm.SPHERE is not None
        self.chemtable = sm.ParserChemTable(self._get_path_table())
        self._init_attrs_from_molecules(path_struct, path_traj)


    # --------------------------------------------------------------------------
    @classmethod
    def from_pqr_data(cls, pqr_data: str):
        with tempfile.NamedTemporaryFile(mode = "w+", suffix = ".pqr", delete = True) as tmp_pqr:
            tmp_pqr.write(pqr_data)
            tmp_pqr.flush()
            return cls(Path(tmp_pqr.name), None)


    # --------------------------------------------------------------------------
    def get_relevant_atoms(self):
        if self.do_ps:
            point = f"{sm.SPHERE.x} {sm.SPHERE.y} {sm.SPHERE.z} {sm.SPHERE.radius}"
            return self.system.select_atoms(
                f"{self.chemtable.selection_query} and point {point}"
            )

        return self.system.select_atoms(self.chemtable.selection_query)


    # --------------------------------------------------------------------------
    def get_relevant_atoms_broad(self, trimming_dist):
        if self.do_ps:
            point = f"{sm.SPHERE.x} {sm.SPHERE.y} {sm.SPHERE.z} {sm.SPHERE.radius + trimming_dist}"
            return self.system.select_atoms(
                f"{self.chemtable.selection_query} and point {point}"
            )

        return self.system.select_atoms(self.chemtable.selection_query)


    # --------------------------------------------------------------------------
    def _infer_box_attributes(self):
        if self.do_ps:
            self.cog = np.array([sm.SPHERE.x, sm.SPHERE.y, sm.SPHERE.z])
            self.min_coords = self.cog - sm.SPHERE.radius
            self.max_coords = self.cog + sm.SPHERE.radius
            self.radius = sm.SPHERE.radius
            return

        super()._infer_box_attributes()


    # --------------------------------------------------------------------------
    def _get_path_table(self) -> Path:
        if sm.PATH_TABLE: return sm.PATH_TABLE

        folder_default_tables = Path("_tables")

        if sm.CURRENT_MOLTYPE == MolType.PROT:
            return vg.resolve_path_package(folder_default_tables / "prot.chem")

        if sm.CURRENT_MOLTYPE == MolType.RNA:
            name = "rna_simple_hb" if sm.DO_SIMPLE_HBONDS_RNA else "rna"
            return vg.resolve_path_package(folder_default_tables / f"{name}.chem")

        raise ValueError(f"No default table for the specified molecular type '{sm.CURRENT_MOLTYPE}'. Please provide a path to a custom table.")


# //////////////////////////////////////////////////////////////////////////////
