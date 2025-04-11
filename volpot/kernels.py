import volpot
import numpy as np

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
        center = np.floor(self.kernel_res / 2) * self.deltas
        coords = np.indices(self.kernel_res).T * self.deltas
        self.shifted_coords = coords - center
        self.dist = volpot.get_norm(self.shifted_coords)

    def link_to_grid(self, grid, grid_origin):
        self.grid = grid
        self.grid_origin = grid_origin
        self.grid_res = np.array(grid.shape)

    def stamp(self, center_stamp_at, multiplication_factor = None):
        if self.grid is None:
            raise ValueError("No grid associated, can't stamp Kernel. Use 'link_to_grid' first.")

        ### infer the position where to stamp the kernel at the big grid
        stamp_orig = center_stamp_at - self.deltas * self.kernel_res / 2
        rel_orig = stamp_orig - self.grid_origin
        idx_orig = np.round(rel_orig / self.deltas).astype(int)
        idx_max = idx_orig + self.kernel_res

        ### skip cases where the kernel would be stamped outside the big grid
        if (idx_max < 0).any(): return
        if (idx_orig > self.grid_res).any(): return

        ### initialize the grid (g_*) and kernel (k_*) indices
        g_i0, g_j0, g_k0 = idx_orig
        g_i1, g_j1, g_k1 = idx_max
        k_i0, k_j0, k_k0 = 0, 0, 0
        k_i1, k_j1, k_k1 = self.kernel_res

        g_rx, g_ry, g_rz = self.grid_res
        k_rx, k_ry, k_rz = self.kernel_res

        ### clamp the indices of both the big grid and the kernel
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


        ### stamp the kernel on the big grid
        if multiplication_factor is None: # multiplication_factor defaults to None (instead of 1) to avoid problems with bool grids
            self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1] += self.kernel[k_i0:k_i1, k_j0:k_j1, k_k0:k_k1]
        else:
            self.grid[g_i0:g_i1, g_j0:g_j1, g_k0:g_k1] += multiplication_factor * self.kernel[k_i0:k_i1, k_j0:k_j1, k_k0:k_k1]


# //////////////////////////////////////////////////////////////////////////////
class SphereKernel(Kernel):
    """For generating simple boolean spheres (e.g. for masks)"""
    def __init__(self, radius, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        self.kernel[self.dist < radius] = 1


# //////////////////////////////////////////////////////////////////////////////
class GaussianKernel(Kernel):
    """For generating univariate gaussian spheres (e.g. for hydrophob)"""
    def __init__(self, mu, sigma, radius, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        self.kernel = volpot.univariate_gaussian(self.dist, mu, sigma)


# //////////////////////////////////////////////////////////////////////////////
class MultiGaussianKernel(Kernel):
    """For generating multivariate gaussian distributions (for hba, hbd, stacking)"""
    def recalculate_kernel(self, normal, mu, cov_inv, isStacking):
        beta_values = volpot.get_angle(self.shifted_coords, normal, isStacking)
        input_mat = np.concatenate(
            (
                np.resize(beta_values, list(beta_values.shape) + [1]),
                np.resize(self.dist,   list(self.dist.shape)   + [1]),
            ),
            axis = 3
        )
        self.kernel = volpot.multivariate_gaussian(input_mat, mu, cov_inv)


# //////////////////////////////////////////////////////////////////////////////
