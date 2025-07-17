import numpy as np
from pathlib import Path
from enum import Enum, auto

import volgrids.vgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class MolType(Enum):
    NONE = auto()
    PROT = auto()
    RNA = auto()
    LIGAND = auto()


# //////////////////////////////////////////////////////////////////////////////
class MolSystemSmiffer(vg.MolSystem):
    def __init__(self, path_struct: Path, path_traj: Path = None):
        self.do_ps = sm.PS_INFO is not None
        self.chemtable = sm.ChemTable(self._get_path_table())
        self._init_attrs_from_molecules(path_struct, path_traj)


    # --------------------------------------------------------------------------
    def get_relevant_atoms(self):
        if self.do_ps:
            radius, xcog, ycog, zcog = sm.PS_INFO
            return self.system.select_atoms(
                f"{self.chemtable.selection_query} and point {xcog} {ycog} {zcog} {radius}"
            )

        return self.system.select_atoms(self.chemtable.selection_query)


    # --------------------------------------------------------------------------
    def get_relevant_atoms_broad(self, trimming_dist):
        if self.do_ps:
            radius, xcog, ycog, zcog = sm.PS_INFO
            return self.system.select_atoms(
                f"{self.chemtable.selection_query} and point {xcog} {ycog} {zcog} {radius + trimming_dist}"
            )

        return self.system.select_atoms(self.chemtable.selection_query)


    # --------------------------------------------------------------------------
    def _infer_box_attributes(self):
        if self.do_ps:
            radius, xcog, ycog, zcog = sm.PS_INFO
            self.cog = np.array([xcog, ycog, zcog])
            self.minCoords = self.cog - radius
            self.maxCoords = self.cog + radius
            self.radius = radius
            return

        super()._infer_box_attributes()


    # --------------------------------------------------------------------------
    def _get_path_table(self) -> Path:
        if sm.PATH_TABLE: return sm.PATH_TABLE

        folder_default_tables = Path("volgrids/smiffer/tables")

        if sm.CURRENT_MOLTYPE == MolType.PROT:
            return vg.resolve_path(folder_default_tables / "prot.chem")

        if sm.CURRENT_MOLTYPE == MolType.RNA:
            return vg.resolve_path(folder_default_tables / "rna.chem")

        raise ValueError(f"No default table for the specified molecular type '{sm.CURRENT_MOLTYPE}'. Please provide a path to a custom table.")


# //////////////////////////////////////////////////////////////////////////////
