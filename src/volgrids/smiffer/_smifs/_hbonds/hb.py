import numpy as np
from abc import ABC, abstractmethod

import volgrids as vg
import volgrids.smiffer as sm

from .utils import get_direction_vector, str_this_residue

# //////////////////////////////////////////////////////////////////////////////
class SmifHBonds(sm.Smif, ABC):
    # --------------------------------------------------------------------------
    @abstractmethod
    def select_antecedent(self, res, res_atoms, hbond_triplet) -> None|np.ndarray:
        return


    # --------------------------------------------------------------------------
    @abstractmethod
    def init_kernel(self):
        return


    # --------------------------------------------------------------------------
    def populate_grid(self):
        self.kernel: vg.Kernel
        self.kernel_args: dict
        self.hbond_getter: callable
        self.init_kernel()

        self.kernel.link_to_grid(self.grid, self.ms.minCoords)
        for pos_antecedent, pos_interactor in self.iter_particles():
            direction = get_direction_vector(vec_origin = pos_antecedent, vec_head = pos_interactor)
            self.kernel.recalculate_kernel(direction, **self.kernel_args)
            self.kernel.stamp(pos_interactor, multiplication_factor = sm.ENERGY_SCALE)


    # --------------------------------------------------------------------------
    def iter_particles(self):
        atoms = self.ms.get_relevant_atoms()
        for res in atoms.residues:
            res_atoms = atoms.select_atoms(str_this_residue(res))
            hbond_triplets = self.hbond_getter(self.ms.chemtable, res.resname)
            if hbond_triplets is None: continue # skip weird residues

            for hbond_triplet in hbond_triplets:
                if not hbond_triplet: continue  # skip residues without HBond pairs

                name_interactor = hbond_triplet[2]

                pos_antecedent = self.select_antecedent(res, res_atoms, hbond_triplet)
                sel_interactor = res_atoms.select_atoms(f"name {name_interactor}")
                if pos_antecedent is None: continue # skip special cases

                if not sel_interactor or len(pos_antecedent) == 0: continue # skip residue artifacts
                pos_interactor = sel_interactor[0].position

                yield pos_antecedent, pos_interactor


# //////////////////////////////////////////////////////////////////////////////
