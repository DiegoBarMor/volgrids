import volgrids as vg
import volgrids.smiffer as sm

from .hb import SmifHBonds
from .triplet import Triplet

# //////////////////////////////////////////////////////////////////////////////
class SmifHBAccepts(SmifHBonds):
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        triplet.set_pos_head(self.res_atoms)

        ############################### TAIL POSITION
        ### special cases for RNA
        if sm.CURRENT_MOLTYPE == sm.MolType.RNA:
            if triplet.interactor == "O3'": # tail points are in different residues
                triplet.set_pos_tail_custom(
                    atoms = self.all_atoms,
                    query_t0 = triplet.str_this_res,
                    query_t1 = triplet.str_next_res
                )
                return

        triplet.set_pos_tail(self.res_atoms)


    # --------------------------------------------------------------------------
    def populate_grid(self):
        self.kernel = vg.KernelGaussianBivariateAngleDist(
            radius = sm.MU_DIST_HBA + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBA,
            deltas = self.ms.deltas, dtype = vg.FLOAT_DTYPE, params = sm.PARAMS_HBA
        )
        self.kernel.link_to_grid(self.grid, self.ms.minCoords)
        self.hbond_getter = sm.ParserChemTable.get_names_hba
        self.process_kernel()


# //////////////////////////////////////////////////////////////////////////////
