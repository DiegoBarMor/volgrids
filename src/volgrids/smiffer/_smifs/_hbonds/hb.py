from abc import ABC, abstractmethod

import volgrids as vg
import volgrids.smiffer as sm

from .triplet import Triplet

# //////////////////////////////////////////////////////////////////////////////
class SmifHBonds(sm.Smif, ABC):
    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel: vg.Kernel = None
        self.kernel_params: vg.ParamsGaussianBivariate = None
        self.hbond_getter: callable


    # --------------------------------------------------------------------------
    @abstractmethod
    def set_triplet_positions(self, triplet: Triplet, all_atoms, res) -> None:
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @abstractmethod
    def init_kernel(self):
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    def populate_grid(self):
        self.init_kernel()
        if self.kernel is None:
            raise ValueError("Kernel must be initialized in the init_kernel() implementation.")
        self.process_kernel()


    # --------------------------------------------------------------------------
    def process_kernel(self):
        self.kernel.link_to_grid(self.grid, self.ms.minCoords)
        for pos_interactor, vec_direction in self.iter_particles():
            self.kernel.recalculate_kernel(
                vec_direction, self.kernel_params, isStacking = False
            )
            self.kernel.stamp(pos_interactor, multiplication_factor = sm.ENERGY_SCALE)


    # --------------------------------------------------------------------------
    def iter_particles(self):
        all_atoms = self.ms.get_relevant_atoms()
        for res in all_atoms.residues:
            triplet_strs = self.hbond_getter(self.ms.chemtable, res.resname)
            if triplet_strs is None: continue # skip weird residues

            for triplet_str in triplet_strs:
                if not triplet_str: continue  # skip residues without HBond pairs

                triplet = Triplet(*triplet_str)
                self.set_triplet_positions(triplet, all_atoms, res)
                vec_direction = triplet.get_direction_vector()

                if (triplet.pos_interactor is None) or (vec_direction is None):
                    continue

                yield triplet.pos_interactor, vec_direction


# //////////////////////////////////////////////////////////////////////////////
