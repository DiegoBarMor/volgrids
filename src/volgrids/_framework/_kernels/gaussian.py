import numpy as np
from abc import ABC

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class KernelGaussian(vg.Kernel, ABC):
    def __init__(self, radius, deltas, dtype, params: "vg.ParamsGaussian"):
        super().__init__(radius, deltas, dtype)
        self.params = params


# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianUnivariateDist(KernelGaussian):
    """For generating univariate gaussian spheres (e.g. for hydrophob)"""
    def __init__(self, radius, deltas, dtype, params: "vg.ParamsGaussianUnivariate"):
        super().__init__(radius, deltas, dtype, params)
        self.kernel = vg.Math.univariate_gaussian(self.dist, params.mu, params.sigma)


# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianBivariateAngleDist(KernelGaussian):
    """For generating multivariate gaussian distributions (for hba, hbd, stacking)"""
    def recalculate_kernel(self, normal, isStacking: bool):
        self.params: "vg.ParamsGaussianBivariate"
        angles = vg.Math.get_angle(
            self.shifted_coords, normal,
            flag_corrections = "stacking" if isStacking else "hbonds"
        )
        input_mat = np.concatenate(
            (
                np.resize(angles,    list(angles.shape)    + [1]),
                np.resize(self.dist, list(self.dist.shape) + [1]),
            ),
            axis = 3
        )
        self.kernel = vg.Math.bivariate_gaussian(input_mat, self.params.mu, self.params.cov_inv)


# //////////////////////////////////////////////////////////////////////////////
