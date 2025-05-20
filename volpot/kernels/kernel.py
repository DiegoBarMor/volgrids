import numpy as np
import volpot as vp

# //////////////////////////////////////////////////////////////////////////////
class Kernel:
    def __init__(self, radius, deltas, dtype):
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
        self.coords = vp.get_coords_array(self.kernel_res, self.deltas)
        self.shifted_coords = self.coords - self.center
        self.dist = vp.get_norm(self.shifted_coords)

    def link_to_grid(self, grid, grid_origin):
        self.grid = grid
        self.grid_origin = grid_origin
        self.grid_res = np.array(grid.shape)

    def stamp(self, center_stamp_at, multiplication_factor = None, operation = "sum"):
        if self.grid is None:
            raise ValueError("No grid associated, can't stamp Kernel. Use 'link_to_grid' first.")

        ##### infer the position where to stamp the kernel at the big grid
        stamp_orig = center_stamp_at - self.deltas * self.kernel_res / 2
        rel_orig = stamp_orig - self.grid_origin
        idx_orig = np.round(rel_orig / self.deltas).astype(int)
        idx_max = idx_orig + self.kernel_res

        ##### skip cases where the kernel would be stamped outside the big grid
        if (idx_max < 0).any(): return
        if (idx_orig > self.grid_res).any(): return

        ##### initialize the grid (g_*) and kernel (k_*) indices
        g_i0, g_j0, g_k0 = idx_orig
        g_i1, g_j1, g_k1 = idx_max
        k_i0, k_j0, k_k0 = 0, 0, 0
        k_i1, k_j1, k_k1 = self.kernel_res

        g_rx, g_ry, g_rz = self.grid_res
        k_rx, k_ry, k_rz = self.kernel_res

        ##### clamp the indices of both the big grid and the kernel
        if g_i0 < 0 and g_i1 >= g_rx: # when the big grid is smaller than the kernel
            k_i0 = -g_i0
            k_i1 = k_i0 + g_rx
            g_i0 = 0
            g_i1 = g_rx
        elif g_i0 < 0:
            g_i0 = 0
            k_i0 = k_rx - g_i1
        elif g_i1 >= g_rx:
            g_i1 = g_rx
            k_i1 = g_rx - g_i0

        if g_j0 < 0 and g_j1 >= g_rx: # when the big grid is smaller than the kernel
            k_j0 = -g_j0
            k_j1 = k_j0 + g_ry
            g_j0 = 0
            g_j1 = g_ry
        elif g_j0 < 0:
            g_j0 = 0
            k_j0 = k_ry - g_j1
        elif g_j1 >= g_ry:
            g_j1 = g_ry
            k_j1 = g_ry - g_j0

        if g_k0 < 0 and g_k1 >= g_rx: # when the big grid is smaller than the kernel
            k_k0 = -g_k0
            k_k1 = k_k0 + g_rz
            g_k0 = 0
            g_k1 = g_rz
        elif g_k0 < 0:
            g_k0 = 0
            k_k0 = k_rz - g_k1
        elif g_k1 >= g_rz:
            g_k1 = g_rz
            k_k1 = g_rz - g_k0


        ##### stamp the kernel on the big grid
        subkernel = self.kernel[k_i0:k_i1, k_j0:k_j1, k_k0:k_k1]

        # multiplication_factor defaults to None (instead of 1) to avoid problems with bool grids
        scaled_subkernel = subkernel if (multiplication_factor is None) else multiplication_factor * subkernel 


        if operation == "sum":
            self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1] += scaled_subkernel

        elif operation == "min":
            self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1] = np.minimum(
                self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1], scaled_subkernel
            )

        elif operation == "max":
            self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1] = np.maximum(
                self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1], scaled_subkernel
            )

        else:
            raise ValueError(f"Unknown operation: {operation}. Use 'sum', 'min' or 'max'.")

# //////////////////////////////////////////////////////////////////////////////
