import numpy as np

import volgrids.vgrids as vg

# ------------------------------------------------------------------------------
def _clamp_indices(g_idx0, g_idx1, k_idx0, k_idx1, g_res, k_res):
    if (g_idx0 < 0) and (g_res <= g_idx1): # kernel would be larger than the big grid
        return 0, g_res, -g_idx0, -g_idx0+g_res

    if g_idx0 < 0: # kernel would start before the big grid
        return 0, g_idx1, k_res-g_idx1, k_idx1

    if g_idx1 >= g_res: # kernel would end after the big grid
        return g_idx0, g_res, k_idx0, g_res-g_idx0

    return g_idx0, g_idx1, k_idx0, k_idx1


# //////////////////////////////////////////////////////////////////////////////
class Kernel:
    def __init__(self, radius, deltas, dtype, operation = "sum"):
        ##### store kernel values
        self.kernel_res = (np.ceil(radius / deltas) * 2 + 1).astype(int)
        self.deltas = deltas
        self.kernel = np.zeros(self.kernel_res, dtype = dtype)

        ##### initialize empty big-grid values (assign them later with link_to_grid)
        self.grid = None
        self.grid_origin = None
        self.grid_res = None

        ##### initizalize auxiliary kernel of distance values
        self.center = np.floor(self.kernel_res / 2) * self.deltas
        self.coords = vg.Math.get_coords_array(self.kernel_res, self.deltas)
        self.shifted_coords = self.coords - self.center
        self.dist = vg.Math.get_norm(self.shifted_coords)

        ##### set operation
        self.operation: callable[np.array, np.array]
        if   operation == "sum": self.operation = np.add
        elif operation == "min": self.operation = np.minimum
        elif operation == "max": self.operation = np.maximum
        else: raise ValueError(f"Unknown operation: {operation}. Use 'sum', 'min' or 'max'.")


    # --------------------------------------------------------------------------
    def link_to_grid(self, grid, grid_origin):
        self.grid = grid
        self.grid_origin = grid_origin
        self.grid_res = np.array(grid.shape)


    # --------------------------------------------------------------------------
    def recalculate_kernel(self, *args, **kwargs):
        return


    # --------------------------------------------------------------------------
    def stamp(self, center_stamp_at, multiplication_factor = None, operation = "sum"):
        if self.grid is None:
            raise ValueError("No grid associated, can't stamp Kernel. Use 'link_to_grid' first.")

        ##### infer the position where to stamp the kernel at the big grid
        stamp_orig = center_stamp_at - self.deltas * self.kernel_res / 2
        rel_orig = stamp_orig - self.grid_origin
        idx_start = np.round(rel_orig / self.deltas).astype(int)
        idx_end = idx_start + self.kernel_res

        ##### skip cases where the kernel would be stamped outside the big grid
        if (idx_end < 0).any(): return
        if (idx_start > self.grid_res).any(): return

        ##### initialize the grid (g_*) and kernel (k_*) indices
        g_i0, g_j0, g_k0 = idx_start
        g_i1, g_j1, g_k1 = idx_end
        k_i0, k_j0, k_k0 = 0, 0, 0
        k_i1, k_j1, k_k1 = self.kernel_res

        g_rx, g_ry, g_rz = self.grid_res
        k_rx, k_ry, k_rz = self.kernel_res

        ##### clamp the indices of both the big grid and the kernel
        g_i0, g_i1, k_i0, k_i1 = _clamp_indices(g_i0, g_i1, k_i0, k_i1, g_rx, k_rx)
        g_j0, g_j1, k_j0, k_j1 = _clamp_indices(g_j0, g_j1, k_j0, k_j1, g_ry, k_ry)
        g_k0, g_k1, k_k0, k_k1 = _clamp_indices(g_k0, g_k1, k_k0, k_k1, g_rz, k_rz)

        ##### stamp the kernel on the big grid
        subkernel = self.kernel[k_i0:k_i1, k_j0:k_j1, k_k0:k_k1]
        subgrid   = self.grid  [g_i0:g_i1, g_j0:g_j1, g_k0:g_k1]

        ### multiplication_factor defaults to None (instead of 1) to avoid problems with bool grids
        scaled_subkernel = subkernel if (multiplication_factor is None) else multiplication_factor * subkernel

        self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1] = self.operation(subgrid, scaled_subkernel)


# //////////////////////////////////////////////////////////////////////////////
