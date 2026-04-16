import numpy as np
import MDAnalysis as mda
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class MolSystem:
    def __init__(self, path_struct: Path, path_traj: Path = None):
        self.molname : str                 # name of the molecule
        self.do_traj : None | bool         # whether this is a trajectory or a single structure (None if no structure is provided)
        self.system  : None | mda.Universe # MDAnalysis Universe object for the molecular system
        self.frame   : None | int          # current frame number (if trajectory is used)
        self.box     : vg.Box

        self.molname = path_struct.stem
        self.do_traj = path_traj is not None

        if self.do_traj:
            self.system = mda.Universe(str(path_struct), str(path_traj))
            self.frame = 0
        else:
            self.system = mda.Universe(str(path_struct))
            self.frame = None

        self.box = self._get_init_box()
        self._enforce_equilateral_grid()


    # --------------------------------------------------------------------------
    def _get_init_box(self) -> vg.Box:
        box = vg.Box(None, None, None, do_init = False)
        box.min_coords = np.min(self.system.coord.positions, axis = 0) - vg.EXTRA_BOX_SIZE
        box.max_coords = np.max(self.system.coord.positions, axis = 0) + vg.EXTRA_BOX_SIZE
        box.infer_deltas_resolution()
        box.infer_radius_and_cog()
        return box



    # --------------------------------------------------------------------------
    def _enforce_equilateral_grid(self):
        if not vg.ENSURE_EQUILATERAL: return

        max_resolution: np.ndarray = np.max(self.box.resolution)
        res_diff = max_resolution - self.box.resolution

        pad_0 = np.ceil (res_diff / 2).astype(int)
        pad_1 = np.floor(res_diff / 2).astype(int)

        self.box.resolution  = np.array([max_resolution, max_resolution, max_resolution], dtype = int)
        self.box.min_coords -= pad_0 * self.box.deltas
        self.box.max_coords += pad_1 * self.box.deltas
        self.box.infer_radius_and_cog()


    # --------------------------------------------------------------------------
    def get_relevant_atoms(self):
        raise NotImplementedError()


# //////////////////////////////////////////////////////////////////////////////
