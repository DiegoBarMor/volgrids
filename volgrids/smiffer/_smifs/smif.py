from abc import abstractmethod
from pathlib import Path
import numpy as np

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Smif(vg.Grid):
    # --------------------------------------------------------------------------
    def __init__(self, ms: "sm.MolSystemSmiffer", init_grid = True, dtype = None):
        super().__init__(ms, init_grid, dtype)
        self.ms: "sm.MolSystemSmiffer" = ms


    # -------------------------------------------------------------------------- INHERITED METHODS
    def __add__(self, other: "Smif|float|int") -> "Smif": return super().__add__(other)
    def __sub__(self, other: "Smif|float|int") -> "Smif": return super().__sub__(other)
    def __abs__(self) -> "Smif": return super().__abs__()

    @classmethod
    def reverse(cls, other: "Smif") -> "Smif": return super().reverse(other)

    def copy(self) -> "Smif":
        obj = self.__class__(self.ms)
        obj.arr = np.copy(self.arr)
        return obj


    # --------------------------------------------------------------------------
    @classmethod
    def from_grid(cls, grid: "vg.Grid", ms: "sm.MolSystemSmiffer") -> "Smif":
        obj = cls(ms, init_grid = False, dtype = grid.dtype)
        obj.arr = np.copy(grid.arr)
        return obj


    # --------------------------------------------------------------------------
    @abstractmethod
    def populate_grid(self):
        raise NotImplementedError("Subclasses of Smif must implement the populate_grid method.")


    # --------------------------------------------------------------------------
    def save_data_smif(self, folder_out: Path, title: str):
        def add_suffix(path: Path, suffix: str) -> Path:
            return Path(str(path) + suffix)

        if self.ms.do_traj:
            path_out = folder_out / f"{self.ms.molname}.{title}.cmap"
            grid_format = vg.GridFormat.CMAP_PACKED # ignores the GRID_FORMAT_OUTPUT config -> CMAP is the only format that supports multiple frames
            cmap_key = f"{self.ms.molname}.{self.ms.frame:04}"

        else:
            path_out = folder_out / f"{self.ms.molname}.{title}"
            grid_format = vg.GridFormat.from_str(sm.GRID_FORMAT_OUTPUT)
            cmap_key = title

            if   grid_format == vg.GridFormat.DX: path_out = add_suffix(path_out, ".dx")
            elif grid_format == vg.GridFormat.MRC: path_out = add_suffix(path_out, ".mrc")
            elif grid_format == vg.GridFormat.CCP4: path_out = add_suffix(path_out, ".ccp4")

            elif grid_format == vg.GridFormat.CMAP:
                path_out = add_suffix(path_out, ".cmap")
                cmap_key = self.ms.molname

            elif grid_format == vg.GridFormat.CMAP_PACKED:
                path_out = folder_out / f"{self.ms.molname}.cmap"
                cmap_key = f"{self.ms.molname}.{title}"

        super().save_data(path_out, grid_format, cmap_key)



# //////////////////////////////////////////////////////////////////////////////
