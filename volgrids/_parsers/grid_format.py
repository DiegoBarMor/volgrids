from enum import Enum, auto

# //////////////////////////////////////////////////////////////////////////////
class GridFormat(Enum):
    DX = auto()
    MRC = auto()
    CCP4 = auto()
    CMAP = auto()
    CMAP_PACKED = auto()

    # --------------------------------------------------------------------------
    @classmethod
    def from_str(cls, s: str) -> "GridFormat":
        s = s.upper()
        if s == "DX": return cls.DX
        if s == "MRC": return cls.MRC
        if s == "CCP4": return cls.CCP4
        if s == "CMAP": return cls.CMAP
        if s == "CMAP_PACKED": return cls.CMAP_PACKED
        raise ValueError(f"Unknown grid format: {s}.")


# //////////////////////////////////////////////////////////////////////////////
