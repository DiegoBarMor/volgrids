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


    # --------------------------------------------------------------------------
    def populate_grid(self):
        self.kernel = vg.KernelGaussianBivariateAngleDist(
            radius = sm.MU_DIST_HBD + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD,
            deltas = self.ms.deltas, dtype = vg.FLOAT_DTYPE, params = sm.PARAMS_HBD
        )
        self.hbond_getter = sm.ChemTable.get_names_hbd
        self.process_kernel()

        self.kernel = vg.KernelGaussianBivariateAngleDist(
            radius = sm.MU_DIST_HBD_FIXED + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD_FIXED,
            deltas = self.ms.deltas, dtype = vg.FLOAT_DTYPE, params = sm.PARAMS_HBD_FIXED
        )
        self.hbond_getter = sm.ChemTable.get_names_hbd_fixed
        self.process_kernel()


# //////////////////////////////////////////////////////////////////////////////
