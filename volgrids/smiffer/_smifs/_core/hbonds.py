from abc import ABC, abstractmethod

import volgrids as vg
import volgrids.smiffer as sm

from .triplet import Triplet
from .smif import Smif

# //////////////////////////////////////////////////////////////////////////////
class SmifHBonds(Smif, ABC):
    def __init__(self, ms: "sm.MolSystem"):
        super().__init__(ms)
        self.kernel: vg.KernelGaussianBivariateAngleDist = None
        self.hbond_getter: callable
        self.all_atoms = self.ms.get_relevant_atoms()
        self.res_atoms = None
        self.processed_interactors = set()


    # --------------------------------------------------------------------------
    @abstractmethod
    def can_be_interactor(self, triplet: Triplet) -> bool:
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @abstractmethod
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid) -> None:
        grid.reset()
        for pos_interactor, vec_direction in self.iter_particles():
            self.kernel.recalculate_kernel(vec_direction, is_stacking = False)
            self.kernel.stamp(grid, pos_interactor, multiply_by = sm.ENERGY_SCALE)


    # --------------------------------------------------------------------------
    def iter_particles(self):
        for triplet in self._iter_triplets():
            self.find_tail_head_positions(triplet)
            vec_direction = triplet.get_direction_vector()

            if (triplet.pos_interactor is None) or (vec_direction is None):
                continue

            yield triplet.pos_interactor, vec_direction


    # --------------------------------------------------------------------------
    def _iter_triplets(self):
        for res in self.all_atoms.residues:
            hbond_tuples = self.hbond_getter(self.ms.chemtable, res.resname)
            if hbond_tuples is None: continue # skip weird residues

            self.processed_interactors.clear()

            for hbond_tuple in hbond_tuples:
                if not hbond_tuple: continue  # skip residues without HBond pairs

                triplet = Triplet(res, *hbond_tuple)
                if not self.can_be_interactor(triplet): continue

                self.res_atoms = self.all_atoms.select_atoms(triplet.str_this_res)
                triplet.set_pos_interactor(self.res_atoms)
                yield triplet


    # ------------------------------------------------------------------------------
    def _has_prev_res(self, triplet: Triplet) -> bool:
        return len(self.all_atoms.select_atoms(triplet.str_prev_res)) > 0


    # ------------------------------------------------------------------------------
    def _has_next_res(self, triplet: Triplet) -> bool:
        return len(self.all_atoms.select_atoms(triplet.str_next_res)) > 0


# //////////////////////////////////////////////////////////////////////////////
