import numpy as np
from abc import ABC

import volgrids as vg
import volgrids.smiffer as sm

from .hb import SmifHBonds
from .triplet import Triplet
from .utils import str_this_residue, str_next_residue

# //////////////////////////////////////////////////////////////////////////////
class SmifHBAccepts(SmifHBonds, ABC):
    def set_triplet_positions(self, triplet: Triplet, all_atoms, res) -> None:
        res_atoms = all_atoms.select_atoms(str_this_residue(res))
        triplet.set_pos_interactor(res_atoms)
        triplet.set_pos_head(res_atoms)

        ############################### TAIL POSITION
        ### special cases for RNA
        if sm.CURRENT_MOLTYPE == sm.MolType.RNA:
            if triplet.interactor_is("O3'"): # tail points are in different residues
                triplet.set_pos_tail_custom(
                    atoms = all_atoms,
                    query_t0 = str_this_residue(res),
                    query_t1 = str_next_residue(res)
                )
                return

        triplet.set_pos_tail(res_atoms)


# //////////////////////////////////////////////////////////////////////////////
class SmifHBARing(SmifHBAccepts):
    def init_kernel(self):
        radius = sm.MU_DIST_HBA + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBA
        self.kernel = vg.KernelGaussianBivariateAngleDist(radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(params = sm.PARAMS_HBA, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hba


# //////////////////////////////////////////////////////////////////////////////
