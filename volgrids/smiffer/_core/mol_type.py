from enum import Enum, auto

# //////////////////////////////////////////////////////////////////////////////
class MolType(Enum):
    NONE = auto()
    PROT = auto()
    RNA = auto()
    MACRO = auto() # PROT + RNA
    LIGAND = auto()

    # --------------------------------------------------------------------------
    @classmethod
    def from_str(cls, s: str):
        s = s.lower()
        if s == "prot"  : return cls.PROT
        if s == "rna"   : return cls.RNA
        if s == "macro" : return cls.MACRO
        if s == "ligand": return cls.LIGAND
        raise ValueError(f"Invalid molecule type string: '{s}'. Expected one of 'prot', 'rna', or 'ligand'.")

    # --------------------------------------------------------------------------
    def is_none(self):   return self == MolType.NONE
    def is_prot(self):   return self == MolType.PROT or self == MolType.MACRO
    def is_rna(self):    return self == MolType.RNA  or self == MolType.MACRO
    def is_ligand(self): return self == MolType.LIGAND


# //////////////////////////////////////////////////////////////////////////////
