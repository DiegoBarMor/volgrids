import numpy as np
from abc import ABC

import volgrids as vg
import volgrids.smiffer as sm

from .hb import SmifHBonds
from .triplet import Triplet
from .utils import str_prev_residue, str_this_residue, str_next_residue

# //////////////////////////////////////////////////////////////////////////////
class SmifHBDonors(SmifHBonds, ABC):
    def set_triplet_positions(self, triplet: Triplet, all_atoms, res) -> None:
        res_atoms = all_atoms.select_atoms(str_this_residue(res))
        triplet.set_pos_interactor(res_atoms)
        triplet.set_pos_head(res_atoms)

        ############################### TAIL POSITION
        ### special cases for protein
        if sm.CURRENT_MOLTYPE == sm.MolType.PROT:
            if triplet.interactor_is("N"): # tail points are in different residues
                triplet.set_pos_tail_custom(
                    atoms = all_atoms,
                    query_t0 = str_prev_residue(res),
                    query_t1 = str_this_residue(res)
                )
                return

        ### special cases for RNA
        if sm.CURRENT_MOLTYPE == sm.MolType.RNA:
            if triplet.interactor_is("O3'"): # donor only if there is no next residue
                sel_next_res = all_atoms.select_atoms(str_next_residue(res))
                if len(sel_next_res) > 0: return

            if triplet.interactor_is("O5'"): # donor only if there is no previous residue
                sel_prev_res = all_atoms.select_atoms(str_prev_residue(res))
                if len(sel_prev_res) > 0: return

        triplet.set_pos_tail(res_atoms)


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
