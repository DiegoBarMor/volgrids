import numpy as np

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianUnivariateDist(vg.Kernel):
    """For generating univariate gaussian spheres (e.g. for hydrophob)"""
    def __init__(self, params: "vg.ParamsGaussianUnivariate", radius, deltas, dtype):
        super().__init__(radius, deltas, dtype)
        self.kernel = vg.Math.univariate_gaussian(self.dist, params.mu, params.sigma)


# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianBivariateAngleDist(vg.Kernel):
    """For generating multivariate gaussian distributions (for hba, hbd, stacking)"""
    def recalculate_kernel(self, normal, params: "vg.ParamsGaussianBivariate", isStacking: bool):
        beta_values = vg.Math.get_angle(
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
        self.kernel = vg.Math.bivariate_gaussian(input_mat, params.mu, params.cov_inv)


# //////////////////////////////////////////////////////////////////////////////
