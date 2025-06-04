import numpy as np
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianUnivariate(vg.Kernel):
    """For generating univariate gaussian spheres (e.g. for hydrophob)"""
    def __init__(self, mu, sigma, radius, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        self.kernel = vg.univariate_gaussian(self.dist, mu, sigma)


# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianMultivariate(vg.Kernel):
    """For generating multivariate gaussian distributions (for hba, hbd, stacking)"""
    def recalculate_kernel(self, normal, mu, cov_inv, isStacking):
        beta_values = vg.get_angle(
            self.shifted_coords, normal,
            flag_corrections = "stacking" if isStacking else "hbonds"
        )
        input_mat = np.concatenate(
            (
                np.resize(beta_values, list(beta_values.shape) + [1]),
                np.resize(self.dist,   list(self.dist.shape)   + [1]),
            ),
            axis = 3
        )
        self.kernel = vg.multivariate_gaussian(input_mat, mu, cov_inv)


# //////////////////////////////////////////////////////////////////////////////
