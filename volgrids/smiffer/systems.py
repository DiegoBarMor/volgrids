import volgrids.vgrids as vg
import volgrids.smiffer as sm
from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////
class SmifferMolecularSystem(vg.MolecularSystem):
    def __init__(self, meta: "vg.ArgsParser"):
        super().__init__(meta)

        folder_tables = Path("volgrids/smiffer/tables")
        path_table = vg.resolve_path(
            folder_tables / f"{'rna' if self.isNucleic else 'prot'}.atoms"
        )
        self.chemtable = sm.ChemTable(path_table)



# //////////////////////////////////////////////////////////////////////////////
