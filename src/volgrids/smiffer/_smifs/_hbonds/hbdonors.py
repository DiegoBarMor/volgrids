import numpy as np
from abc import ABC

import volgrids as vg
import volgrids.smiffer as sm

from .hb import SmifHBonds
from .utils import safe_return_coords, \
    str_prev_residue, str_this_residue, str_next_residue

# //////////////////////////////////////////////////////////////////////////////
class SmifHBDonors(SmifHBonds, ABC):
    def select_antecedent(self, res, res_atoms, hbond_triplet) -> None|np.ndarray:
        name_antecedent_0, name_antecedent_1, name_interactor, do_infer_H = hbond_triplet
        atoms = self.ms.get_relevant_atoms()

        ##### pseudo-antecedents
        if name_antecedent_1:
            ### special case for terminal aminoacids
            if hbond_triplet == ("C", "CA", "N", False):
                return safe_return_coords(atoms,
                    f"({str_prev_residue(res)} and name C ) or " +\
                    f"({str_this_residue(res)} and name CA)"
                )

            if do_infer_H:
                atom_ref        = res_atoms.select_atoms(f"name {name_antecedent_0}")
                atom_antecedent = res_atoms.select_atoms(f"name {name_antecedent_1}")
                atom_hbond      = res_atoms.select_atoms(f"name {name_interactor  }")
                if len(atom_ref) != 1 or len(atom_antecedent) != 1 or len(atom_hbond) != 1:
                    return

                direction = atom_antecedent[0].position - atom_ref[0].position
                return atom_hbond[0].position - direction

            return safe_return_coords(res_atoms, f"name {name_antecedent_0} {name_antecedent_1}")

        ##### standard antecedents
        ## special case for RNA, it's a donor only if there is no next residue
        if name_interactor == "O3'":
            sel_next_res = atoms.select_atoms(str_next_residue(res))
            if len(sel_next_res) > 0: return

        ## special case for RNA, it's a donor only if there is no previous residue
        if name_interactor == "O5'":
            sel_prev_res = atoms.select_atoms(str_prev_residue(res))
            if len(sel_prev_res) > 0: return

        ## other standard antecedent cases
        return safe_return_coords(res_atoms, f"name {name_antecedent_0}")


# //////////////////////////////////////////////////////////////////////////////
class SmifHBDRing(SmifHBDonors):
    def init_kernel(self):
        radius = sm.MU_HBD[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD
        self.kernel = vg.KernelGaussianMultivariate(radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBD, cov_inv = sm.COV_INV_HBD, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd


# //////////////////////////////////////////////////////////////////////////////
class SmifHBDCone(SmifHBDonors):
    def init_kernel(self):
        radius = sm.MU_HBD_FIXED[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD_FIXED
        self.kernel = vg.KernelGaussianMultivariate(radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBD_FIXED, cov_inv = sm.COV_INV_HBD_FIXED, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd_fixed


# //////////////////////////////////////////////////////////////////////////////
