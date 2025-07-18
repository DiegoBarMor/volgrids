import numpy as np
from abc import ABC

import volgrids as vg
import volgrids.smiffer as sm

from .hb import SmifHBonds
from .utils import safe_return_coords, \
    str_this_residue, str_next_residue

# //////////////////////////////////////////////////////////////////////////////
class SmifHBAccepts(SmifHBonds, ABC):
    def select_antecedent(self, res, res_atoms, hbond_triplet) -> None|np.ndarray:
        name_antecedent_0, name_antecedent_1, name_interactor, _ = hbond_triplet
        atoms = self.ms.get_relevant_atoms()

        ##### standard antecedents
        if not name_antecedent_1:
            return safe_return_coords(res_atoms, f"name {name_antecedent_0}")

        ##### pseudo-antecedents
        ### special case for RNA, needs to check next residue
        if name_interactor == "O3'":
            return safe_return_coords(atoms,
                f"({str_this_residue(res)} and name {name_antecedent_0}) or" +\
                f"({str_next_residue(res)} and name {name_antecedent_1})"
            )

        ### other pseudo-antecedent cases
        return safe_return_coords(res_atoms, f"name {name_antecedent_0} {name_antecedent_1}")


# //////////////////////////////////////////////////////////////////////////////
class SmifHBARing(SmifHBAccepts):
    def init_kernel(self):
        radius = sm.MU_HBA[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBA
        self.kernel = vg.KernelGaussianMultivariate(radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBA, cov_inv = sm.COV_INV_HBA, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hba


# //////////////////////////////////////////////////////////////////////////////
