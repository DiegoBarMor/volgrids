import numpy as np

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Cavities:
    # --------------------------------------------------------------------------
    @staticmethod
    def find_cavities_naive(occupancy: vg.Grid) -> vg.Grid:
        arr: np.ndarray = occupancy.grid.astype(bool)

        xsurf = np.zeros_like(arr, dtype = bool)
        ysurf = np.zeros_like(arr, dtype = bool)
        zsurf = np.zeros_like(arr, dtype = bool)
        xsurf[1:,:,:] = arr[1:,:,:] ^ arr[:-1,:,:]
        ysurf[:,1:,:] = arr[:,1:,:] ^ arr[:,:-1,:]
        zsurf[:,:,1:] = arr[:,:,1:] ^ arr[:,:,:-1]

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

        available_x0 = (xrange                  < xsurf_start )
        available_x1 = (xrange                  >= xsurf_end  )
        available_y0 = (yrange.transpose(1,0,2) <  ysurf_start).transpose(1,0,2)
        available_y1 = (yrange.transpose(1,0,2) >= ysurf_end  ).transpose(1,0,2)
        available_z0 = (zrange.transpose(2,0,1) <  zsurf_start).transpose(1,2,0)
        available_z1 = (zrange.transpose(2,0,1) >= zsurf_end  ).transpose(1,2,0)

        available_x = (available_x0 | available_x1).astype(int)
        available_y = (available_y0 | available_y1).astype(int)
        available_z = (available_z0 | available_z1).astype(int)

        cavities = 3 - (available_x + available_y + available_z)
        cavities[arr] = 0

        vgrid = vg.Grid(occupancy.ms, init_grid = False)
        vgrid.grid = np.copy(cavities)
        return vgrid


# //////////////////////////////////////////////////////////////////////////////
