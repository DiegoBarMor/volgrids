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
class SmifferMolecularSystem(vg.MolecularSystem):
    def __init__(self, meta: "sm.SmifferArgsParser"):
        super().__init__(meta)

        self.chemtable = sm.ChemTable(self._get_path_table(meta))

        if meta.do_ps:
            radius,xcog,ycog,zcog = meta.ps_info
            self.relevant_atoms = self.system.select_atoms(
                f"{self.chemtable.selection_query} and point {xcog} {ycog} {zcog} {radius}"
            )
        else:
            self.relevant_atoms = self.system.select_atoms(self.chemtable.selection_query)


    # --------------------------------------------------------------------------
    def _get_path_table(self, meta: "sm.SmifferArgsParser") -> Path:
        if meta.path_table: return meta.path_table

        folder_default_tables = Path("smiffer/tables")

        if meta.moltype == MolType.PROT:
            return vg.resolve_path(folder_default_tables / "prot.atoms")

        if meta.moltype == MolType.RNA:
            return vg.resolve_path(folder_default_tables / "rna.atoms")

        raise ValueError(f"No default table for the specified molecular type '{meta.moltype}'. Please provide a path to a custom table.")


# //////////////////////////////////////////////////////////////////////////////
