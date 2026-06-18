from abc import abstractmethod
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Smif:
    # --------------------------------------------------------------------------
    def __init__(self, ms: "sm.MolSystem"):
        self.ms: "sm.MolSystem" = ms


    # --------------------------------------------------------------------------
    @abstractmethod
    def populate_grid(self, grid: vg.Grid) -> None:
        """
        Inherited versions of `populate_grid` should start by calling `grid.reset()`.
        `grid.dirty = True` will be automatically set by the calls to `Kernel.stamp()`.
        """
        grid.reset()
        raise NotImplementedError("Subclasses of Smif must implement the populate_grid method.")


    # --------------------------------------------------------------------------
    @staticmethod
    def save_data(grid: vg.Grid, ms: sm.MolSystem, folder_out: Path, title: str):
        def add_suffix(path: Path, suffix: str) -> Path:
            return Path(str(path) + suffix)

        if ms.do_traj:
            path_out = folder_out / f"{ms.molname}.{title}.cmap"
            grid_format = vg.GridFormat.CMAP_PACKED # ignores the GRID_FORMAT_OUTPUT config -> CMAP is the only format that supports multiple frames
            key_cmap = f"{ms.molname}.{ms.frame:04}"

        else:
            path_out = folder_out / f"{ms.molname}.{title}"
            grid_format = vg.GridFormat.from_str(sm.GRID_FORMAT_OUTPUT)
            key_cmap = title

            if   grid_format == vg.GridFormat.DX: path_out = add_suffix(path_out, ".dx")
            elif grid_format == vg.GridFormat.BIN: path_out = add_suffix(path_out, ".bin")
            elif grid_format == vg.GridFormat.MRC: path_out = add_suffix(path_out, ".mrc")
            elif grid_format == vg.GridFormat.CCP4: path_out = add_suffix(path_out, ".ccp4")

            elif grid_format == vg.GridFormat.CMAP:
                path_out = add_suffix(path_out, ".cmap")
                key_cmap = ms.molname

            elif grid_format == vg.GridFormat.CMAP_PACKED:
                path_out = folder_out / f"{ms.molname}.cmap"
                key_cmap = f"{ms.molname}.{title}"

        grid.save(path_out, grid_format, key_cmap)



# //////////////////////////////////////////////////////////////////////////////
