import numpy as np
from abc import ABC

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class KernelGaussian(vg.Kernel, ABC):
    def __init__(self,
        radius, deltas, dtype, params: "vg.ParamsGaussian",
        kop: vg.KOperation = vg.KOperation.ADD
    ):
        super().__init__(radius, deltas, dtype, kop)
        self.params = params


# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianUnivariateDist(KernelGaussian):
    """For generating univariate gaussian spheres (e.g. for hydrophob)"""
    def __init__(self,
        radius, deltas, dtype, params: "vg.ParamsGaussianUnivariate",
        kop: vg.KOperation = vg.KOperation.ADD
    ):
        super().__init__(radius, deltas, dtype, params, kop)
        self.arr = vg.Math.univariate_gaussian(self.dists, params.mu, params.sigma)


# //////////////////////////////////////////////////////////////////////////////
class KernelGaussianBivariateAngleDist(KernelGaussian):
    """For generating multivariate gaussian distributions (for hba, hbd, stacking)"""
    def recalculate_kernel(self, normal, is_stacking: bool):
        self.params: "vg.ParamsGaussianBivariate"
        input_mat = self.get_input_mat(self.coords, self.dists, normal, is_stacking)
        self.arr = vg.Math.bivariate_gaussian(input_mat, self.params.mu, self.params.cov_inv)


    # --------------------------------------------------------------------------
    @staticmethod
    def get_input_mat(
        coords: np.ndarray, dists: np.ndarray, normal: np.ndarray, is_stacking: bool
    ) -> np.ndarray:
        angles = vg.Math.get_angle(
            coords, normal,
            flag_corrections = "stacking" if is_stacking else "hbonds"
        )
        return np.concatenate(
            (
                np.resize(angles, list(angles.shape) + [1]),
                np.resize(dists,  list(dists.shape)  + [1]),
            ),
            axis = 3
        )


# //////////////////////////////////////////////////////////////////////////////
