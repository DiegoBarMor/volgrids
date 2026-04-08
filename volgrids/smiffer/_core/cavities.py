import numpy as np

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Cavities:
    # --------------------------------------------------------------------------
    @staticmethod
    def find_cavities_naive(occupancy: vg.Grid) -> vg.Grid:
        arr: np.ndarray = occupancy.grid.astype(bool)

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
        cavities = 3 - (available_x + available_y + available_z)
        cavities[arr] = 0

        vgrid = vg.Grid(occupancy.ms, init_grid = False)
        vgrid.grid = np.copy(cavities)
        return vgrid

    # --------------------------------------------------------------------------
    @staticmethod
    def find_cavities_naive_double_pass(occupancy: vg.Grid) -> vg.Grid:
        ### Step 1: consider the original occupancy grid, and a copy of it rotated by 45 degrees in all dimensions
        occupy_00deg = occupancy
        occupy_45deg = vg.Grid(occupancy.ms, init_grid = False)
        occupy_45deg.grid = vg.Math.rotate_3d(occupancy.grid, 45, 45, 45)

        ### Step 2: find the cavities in both grids with the naive method.
        cavities_00deg = Cavities.find_cavities_naive(occupy_00deg)
        cavities_45deg = Cavities.find_cavities_naive(occupy_45deg)

        ### Step 3: revert the rotated cavities grid back to the original orientation
        cavities_45deg.grid = vg.Math.rotate_3d(cavities_45deg.grid, 45, 45, 45, reverse = True)

        ### Step 4: average the two cavity grids. Discrete values will now be 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0.
        cavities_00deg.grid = (cavities_00deg.grid + cavities_45deg.grid) / 2
        return cavities_00deg


# //////////////////////////////////////////////////////////////////////////////
