import numpy as np

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class CavityFinder:
    CAVITIES_MAX_VAL = 3

    # --------------------------------------------------------------------------
    def __init__(self):
        self.grid: vg.Grid = None


    # --------------------------------------------------------------------------
    def has_data(self) -> bool:
        return self.grid is not None


    # --------------------------------------------------------------------------
    def populate_cavities_grid(self, occupancy: vg.Grid):
        func: callable = self._find_cavities_naive_multi_pass \
            if sm.CAVITIES_NPASSES > 1 else self._find_cavities_naive
        self.grid = func(occupancy)


    # --------------------------------------------------------------------------
    def apply_cavities_weighting(self, smif: "sm.Smif"):
        if sm.CAVITIES_WEIGHT == 0.0: return
        if not self.has_data(): return
        if not self.grid.has_equivalent_box(smif.grid.box):
            print("WARNING: Cavity grid and smif grid do not have the same box. Cavity weighting aborted.")
            return

        smif.grid.arr *= (1 + self.grid.arr * sm.CAVITIES_WEIGHT)


    # --------------------------------------------------------------------------
    @classmethod
    def _find_cavities_naive(cls, occupancy: vg.Grid) -> vg.Grid:
        arr: np.ndarray = occupancy.arr.astype(bool)

        ### STEP 1: find the occupancy surface with XOR
        xsurf = np.zeros_like(arr, dtype = bool)
        ysurf = np.zeros_like(arr, dtype = bool)
        zsurf = np.zeros_like(arr, dtype = bool)
        xsurf[1:,:,:] = arr[1:,:,:] ^ arr[:-1,:,:]
        ysurf[:,1:,:] = arr[:,1:,:] ^ arr[:,:-1,:]
        zsurf[:,:,1:] = arr[:,:,1:] ^ arr[:,:,:-1]

        ### STEP 2: find the indices for the beginning and end of the surface along every dimension
        xrange = np.broadcast_to(np.arange(arr.shape[0])[:,None,None], arr.shape)
        yrange = np.broadcast_to(np.arange(arr.shape[1])[None,:,None], arr.shape)
        zrange = np.broadcast_to(np.arange(arr.shape[2])[None,None,:], arr.shape)

        xsurf_start = np.argmax(xsurf, axis = 0)
        ysurf_start = np.argmax(ysurf, axis = 1)
        zsurf_start = np.argmax(zsurf, axis = 2)
        xsurf_end   = xrange.shape[0] - np.argmax(xsurf[::-1,:,:], axis = 0) - 1
        ysurf_end   = yrange.shape[1] - np.argmax(ysurf[:,::-1,:], axis = 1) - 1
        zsurf_end   = zrange.shape[2] - np.argmax(zsurf[:,:,::-1], axis = 2) - 1

        xsurf_start[xsurf_start == 0] = arr.shape[0] # 0 implies no surface, so set the volume start to end of the dimension
        ysurf_start[ysurf_start == 0] = arr.shape[1]
        zsurf_start[zsurf_start == 0] = arr.shape[2]

        ### STEP 3: populate with 1 the volume before and after the surface limits. Repeat for every dimension
        available_x0 = (xrange                  <  xsurf_start)
        available_x1 = (xrange                  >= xsurf_end  )
        available_y0 = (yrange.transpose(1,0,2) <  ysurf_start).transpose(1,0,2)
        available_y1 = (yrange.transpose(1,0,2) >= ysurf_end  ).transpose(1,0,2)
        available_z0 = (zrange.transpose(2,0,1) <  zsurf_start).transpose(1,2,0)
        available_z1 = (zrange.transpose(2,0,1) >= zsurf_end  ).transpose(1,2,0)

        available_x = (available_x0 | available_x1).astype(int)
        available_y = (available_y0 | available_y1).astype(int)
        available_z = (available_z0 | available_z1).astype(int)

        ### STEP 4: sum 3 dimesions and invert the values, so that deeper cavities have higher values.
        ### Set the occupied volume to 0. At the end, this grid has discrete values 0,1,2,3.
        cavities = cls.CAVITIES_MAX_VAL - (available_x + available_y + available_z)
        cavities[arr] = 0

        vgrid = vg.Grid(occupancy.box, init_grid = False)
        vgrid.arr = np.copy(cavities)
        return vgrid


    # --------------------------------------------------------------------------
    @staticmethod
    def _find_cavities_naive_multi_pass(occupancy: vg.Grid) -> vg.Grid:
        def get_cavities_rot(angle: float):
            occupy_rot = vg.Grid(occupancy.box, init_grid = False)
            occupy_rot.arr = vg.Math.rotate_3d(occupancy.arr, angle, angle, angle)
            cavities_rot = CavityFinder._find_cavities_naive(occupy_rot)
            cavities_rot.arr = vg.Math.rotate_3d(cavities_rot.arr, angle, angle, angle, reverse = True)
            return cavities_rot

        cavities_orig = CavityFinder._find_cavities_naive(occupancy)

        cavities_rot = (
            get_cavities_rot(i*90 / sm.CAVITIES_NPASSES)
            for i in range(1, sm.CAVITIES_NPASSES)
        )

        cavities_orig.arr = (cavities_orig.arr + sum(cav.arr for cav in cavities_rot)) / sm.CAVITIES_NPASSES
        return cavities_orig


# //////////////////////////////////////////////////////////////////////////////
