from abc import ABC, abstractmethod

import volgrids as vg
import volgrids.smiffer as smf

from .triplet import Triplet
from .smif import Smif

# //////////////////////////////////////////////////////////////////////////////
class SmifHBonds(Smif, ABC):
    def __init__(self, mm: "smf.MoleculeManager"):
        super().__init__(mm)
        self.kernel: vg.KernelGaussianBivariateAngleDist = None
        self.dict_triplets: dict[str, list[Triplet]] = {}
        self.atoms = self.mm.get_atoms_insphere()
        self.chains = self.atoms.split_chains()
        self.processed_interactors: set[str] = set()


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
            self.kernel.stamp(grid, pos_interactor, multiply_by = vg.CFG.param_hb_scale)


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
        for chain in self.chains:
            residues = chain.split_residues()

            for i,res in enumerate(residues):
                resname = res[0].resname
                triplets: list[Triplet]|None = self.dict_triplets.get(resname)
                if triplets is None: continue # skip weird residues

                self.processed_interactors.clear()

                for triplet in triplets:
                    res_prev = residues[i-1] if i > 0 else None
                    res_next = residues[i+1] if i < len(residues)-1 else None
                    triplet.resname = resname
                    triplet.residue_prev = res_prev
                    triplet.residue_this = res
                    triplet.residue_next = res_next

                    if not self.can_be_interactor(triplet): continue

                    triplet.reset_posits()
                    triplet.set_pos_interactor()
                    yield triplet


# //////////////////////////////////////////////////////////////////////////////
