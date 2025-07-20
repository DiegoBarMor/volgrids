from abc import ABC, abstractmethod

import volgrids as vg
import volgrids.smiffer as sm

from .triplet import Triplet

# //////////////////////////////////////////////////////////////////////////////
class SmifHBonds(sm.Smif, ABC):
    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel: vg.KernelGaussianBivariateAngleDist = None
        self.hbond_getter: callable
        self.all_atoms = self.ms.get_relevant_atoms()
        self.res_atoms = None
        self.processed_interactors = set()

    # --------------------------------------------------------------------------
    @abstractmethod
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    def process_kernel(self):
        for pos_interactor, vec_direction in self.iter_particles():
            self.kernel.recalculate_kernel(vec_direction, isStacking = False)
            self.kernel.stamp(pos_interactor, multiplication_factor = sm.ENERGY_SCALE)


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
            triplet_strs = self.hbond_getter(self.ms.chemtable, res.resname)
            if triplet_strs is None: continue # skip weird residues

            self.processed_interactors.clear()

            for triplet_str in triplet_strs:
                if not triplet_str: continue  # skip residues without HBond pairs

                triplet = Triplet(res, *triplet_str)
                self.res_atoms = self.all_atoms.select_atoms(triplet.str_this_res)
                triplet.set_pos_interactor(self.res_atoms)
                yield triplet


# //////////////////////////////////////////////////////////////////////////////
