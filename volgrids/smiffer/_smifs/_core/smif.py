from abc import abstractmethod
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Smif:
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
    def save_data(grid: vg.Grid, ms: sm.MolSystem, path_out: Path, key_out: str) -> None:
        # trajectory ignores the GRID_FORMAT_OUTPUT config -> CMAP is the only format that supports multiple frames
        fmt = vg.GridFormat.CMAP if ms.enforce_cmap_output \
            else vg.GridFormat.from_str(sm.GRID_FORMAT_OUTPUT)
        if ms.do_traj: key_out += f".{ms.frame:04}"
        grid.save(path_out, fmt, key_out)


# //////////////////////////////////////////////////////////////////////////////
