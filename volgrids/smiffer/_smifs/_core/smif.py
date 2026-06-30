from abc import abstractmethod
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as smf

# //////////////////////////////////////////////////////////////////////////////
class Smif:
    def __init__(self, mm: "smf.MoleculeManager"):
        self.mm: "smf.MoleculeManager" = mm


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
    def save_data(grid: vg.Grid, mm: smf.MoleculeManager, path_out: Path, key_out: str) -> None:
        # trajectory ignores the OUT_FORMAT config -> CMAP is the only format that supports multiple frames
        fmt = vg.GridFormat.CMAP if mm.enforce_cmap_output \
            else vg.GridFormat.from_str(vg.CFG.out_format)
        if mm.do_traj: key_out += f".{mm.frame:04}"
        grid.save(path_out, fmt, key_out)


# //////////////////////////////////////////////////////////////////////////////
