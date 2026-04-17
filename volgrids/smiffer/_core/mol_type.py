from enum import Enum, auto

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
