from abc import abstractmethod
from pathlib import Path
import numpy as np

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Smif:
    # --------------------------------------------------------------------------
    def __init__(self, ms: "sm.MolSystemSmiffer"):
        self.grid = vg.Grid(ms.box)
        self.ms: "sm.MolSystemSmiffer" = ms


    # --------------------------------------------------------------------------
    @abstractmethod
    def populate_grid(self):
        raise NotImplementedError("Subclasses of Smif must implement the populate_grid method.")


    # --------------------------------------------------------------------------
    @staticmethod
    def save_data_smif(grid: vg.Grid, ms: sm.MolSystemSmiffer, folder_out: Path, title: str):
        def add_suffix(path: Path, suffix: str) -> Path:
            return Path(str(path) + suffix)

        if ms.do_traj:
            path_out = folder_out / f"{ms.molname}.{title}.cmap"
            grid_format = vg.GridFormat.CMAP_PACKED # ignores the GRID_FORMAT_OUTPUT config -> CMAP is the only format that supports multiple frames
            cmap_key = f"{ms.molname}.{ms.frame:04}"

        else:
            path_out = folder_out / f"{ms.molname}.{title}"
            grid_format = vg.GridFormat.from_str(sm.GRID_FORMAT_OUTPUT)
            cmap_key = title

            if   grid_format == vg.GridFormat.DX: path_out = add_suffix(path_out, ".dx")
            elif grid_format == vg.GridFormat.MRC: path_out = add_suffix(path_out, ".mrc")
            elif grid_format == vg.GridFormat.CCP4: path_out = add_suffix(path_out, ".ccp4")

            elif grid_format == vg.GridFormat.CMAP:
                path_out = add_suffix(path_out, ".cmap")
                cmap_key = ms.molname

            elif grid_format == vg.GridFormat.CMAP_PACKED:
                path_out = folder_out / f"{ms.molname}.cmap"
                cmap_key = f"{ms.molname}.{title}"

        grid.save_data(path_out, grid_format, cmap_key)



# //////////////////////////////////////////////////////////////////////////////
