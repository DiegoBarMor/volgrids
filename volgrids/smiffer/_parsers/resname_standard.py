# //////////////////////////////////////////////////////////////////////////////
class ResnameStandard:
    ### [WIP] add a more comprehensive mapping of aliases
    ### see: https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/core/selection.py
    ALIASES_RNA = {
         "U": "U",  "C": "C",  "A": "A",  "G": "G", # standard RNA residue names
        "RU": "U", "RC": "C", "RA": "A", "RG": "G", # alternative RNA residue names
    }
    ALIASES_DNA = {
        "DT": "DT", "DC": "DC", "DA": "DA", "DG": "DG", # standard DNA residue names
    }
    ALIASES_PROT = {
        "ALA": "ALA", "ARG": "ARG", "ASN": "ASN", "ASP": "ASP", "CYS": "CYS",
        "GLU": "GLU", "GLN": "GLN", "GLY": "GLY", "HIS": "HIS", "ILE": "ILE",
        "LEU": "LEU", "LYS": "LYS", "MET": "MET", "PHE": "PHE", "PRO": "PRO",
        "SER": "SER", "THR": "THR", "TRP": "TRP", "TYR": "TYR", "VAL": "VAL",
    }

    # --------------------------------------------------------------------------
    @classmethod
    def standardize(cls, resname: str) -> str:
        #### [TODO] could eventually run the standardization once as an initial preprocessing of the PDB contents,
        #### perhaps when replacing mda as main PDB parser...
        if cls.is_rna (resname): return cls.ALIASES_RNA [resname]
        if cls.is_dna (resname): return cls.ALIASES_DNA [resname]
        if cls.is_prot(resname): return cls.ALIASES_PROT[resname]
        return resname

    # --------------------------------------------------------------------------
    @classmethod
    def is_rna(cls, resname: str) -> bool:
        return resname in cls.ALIASES_RNA

    # --------------------------------------------------------------------------
    @classmethod
    def is_dna(cls, resname: str) -> bool:
        return resname in cls.ALIASES_DNA

    # --------------------------------------------------------------------------
    @classmethod
    def is_nucleic(cls, resname: str) -> bool:
        return cls.is_rna(resname) or cls.is_dna(resname)

    # --------------------------------------------------------------------------
    @classmethod
    def is_prot(cls, resname: str) -> bool:
        return resname in cls.ALIASES_PROT


# //////////////////////////////////////////////////////////////////////////////
