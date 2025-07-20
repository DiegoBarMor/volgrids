import numpy as np
from dataclasses import dataclass, field

# //////////////////////////////////////////////////////////////////////////////
@dataclass
class ParamsGaussianUnivariate:
    mu: float
    sigma: float


# //////////////////////////////////////////////////////////////////////////////
@dataclass
class ParamsGaussianBivariate:
    mu_0: float
    mu_1: float
    cov_00: float
    cov_01: float
    cov_10: float
    cov_11: float
    mu: np.ndarray = field(init = False)
    cov: np.ndarray = field(init = False)
    cov_inv: np.ndarray = field(init = False)

    # --------------------------------------------------------------------------
    def __post_init__(self):
        self.mu = np.array([self.mu_0, self.mu_1])
        self.cov = np.array([
            [self.cov_00, self.cov_01],
            [self.cov_10, self.cov_11],
        ])
        self.cov_inv = np.linalg.inv(self.cov)


# //////////////////////////////////////////////////////////////////////////////
