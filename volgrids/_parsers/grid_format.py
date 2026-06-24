from enum import Enum, auto

# //////////////////////////////////////////////////////////////////////////////
class GridFormat(Enum):
    DX = auto()
    BIN = auto()
    MRC = auto()
    CCP4 = auto()
    CMAP = auto()

    # --------------------------------------------------------------------------
    @classmethod
    def from_str(cls, s: str) -> "GridFormat":
        s = s.upper()
        if s == "DX": return cls.DX
        if s == "BIN": return cls.BIN
        if s == "MRC": return cls.MRC
        if s == "CCP4": return cls.CCP4
        if s == "CMAP": return cls.CMAP
        known_formats = ", ".join(fmt.name for fmt in cls)
        raise ValueError(f"Unknown grid format: {s}. Known formats: {known_formats}.")


    # --------------------------------------------------------------------------
    def suffix(self) -> str:
        return self.name.lower()


# //////////////////////////////////////////////////////////////////////////////
