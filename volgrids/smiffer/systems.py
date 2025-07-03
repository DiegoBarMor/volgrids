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


# //////////////////////////////////////////////////////////////////////////////
class SmifferMolecularSystem(vg.MolecularSystem):
    def __init__(self, path_structure: Path, path_trajectory: Path):
        self.do_ps = sm.PS_INFO is not None
        self.chemtable = sm.ChemTable(self._get_path_table())

        super().__init__(path_structure, path_trajectory)

        if self.do_ps:
            radius,xcog,ycog,zcog = sm.PS_INFO
            self.relevant_atoms = self.system.select_atoms(
                f"{self.chemtable.selection_query} and point {xcog} {ycog} {zcog} {radius}"
            )
        else:
            self.relevant_atoms = self.system.select_atoms(self.chemtable.selection_query)


    # --------------------------------------------------------------------------
    def get_relevant_atoms_broad(self, trimming_dist):
        if not self.do_ps:
            return self.relevant_atoms

        xcog, ycog, zcog = self.cog
        return self.system.select_atoms(
            f"{self.chemtable.selection_query} and point {xcog} {ycog} {zcog} {self.radius + trimming_dist}"
        )


    # --------------------------------------------------------------------------
    def _init_box_attributes(self):
        if self.do_ps:
            radius,xcog,ycog,zcog = sm.PS_INFO
            self.cog = np.array([xcog, ycog, zcog])
            self.minCoords = self.cog - radius
            self.maxCoords = self.cog + radius
            self.radius = radius

        else:
            super()._init_box_attributes()


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
