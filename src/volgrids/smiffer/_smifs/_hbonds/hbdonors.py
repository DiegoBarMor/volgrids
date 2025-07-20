from abc import ABC

import volgrids as vg
import volgrids.smiffer as sm

from .hb import SmifHBonds
from .triplet import Triplet

# ------------------------------------------------------------------------------
def _has_prev_res(atoms, triplet: Triplet) -> bool:
    return len(atoms.select_atoms(triplet.str_prev_res)) > 0

# ------------------------------------------------------------------------------
def _has_next_res(atoms, triplet: Triplet) -> bool:
    return len(atoms.select_atoms(triplet.str_next_res)) > 0


# //////////////////////////////////////////////////////////////////////////////
class SmifHBDonors(SmifHBonds, ABC):
    # --------------------------------------------------------------------------
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        if triplet.pos_head is not None: # head position is already set for succesful sm.USE_STRUCTURE_HYDROGENS iterations
            return

        triplet.set_pos_head(self.res_atoms)

        ############################### TAIL POSITION
        ### special cases for protein
        if sm.CURRENT_MOLTYPE == sm.MolType.PROT:
            if triplet.resname == "PRO": # donor only if there is no previous residue
                if _has_prev_res(self.all_atoms, triplet): return

            elif triplet.interactor == "N": # tail points are in different residues
                if _has_prev_res(self.all_atoms, triplet):
                    triplet.set_pos_tail_custom( # N of peptide bond
                        atoms = self.all_atoms,
                        query_t0 = triplet.str_prev_res,
                        query_t1 = triplet.str_this_res
                    )
                    self.kernel = self._kernel_alt
                    return

                triplet.set_pos_tail_custom( # N of N-terminus
                    atoms = self.all_atoms,
                    query_t0 = f"{triplet.str_this_res} and name CA",
                    query_t1 = f"{triplet.str_this_res} and name CA"
                )
                self.kernel = self._kernel_std
                return


        ### special cases for RNA
        if sm.CURRENT_MOLTYPE == sm.MolType.RNA:
            if triplet.interactor == "O3'": # donor only if there is no next residue
                if _has_next_res(self.all_atoms, triplet): return

            elif triplet.interactor == "O5'": # donor only if there is no previous residue
                if _has_prev_res(self.all_atoms, triplet): return

        triplet.set_pos_tail(self.res_atoms)


    # --------------------------------------------------------------------------
    def populate_grid(self):
        _kernel_hbd = vg.KernelGaussianBivariateAngleDist(
            radius = sm.MU_DIST_HBD + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD,
            deltas = self.ms.deltas, dtype = vg.FLOAT_DTYPE, params = sm.PARAMS_HBD
        )
        _kernel_hbd.link_to_grid(self.grid, self.ms.minCoords)

        _kernel_hbd_fixed = vg.KernelGaussianBivariateAngleDist(
            radius = sm.MU_DIST_HBD_FIXED + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD_FIXED,
            deltas = self.ms.deltas, dtype = vg.FLOAT_DTYPE, params = sm.PARAMS_HBD_FIXED
        )
        _kernel_hbd_fixed.link_to_grid(self.grid, self.ms.minCoords)

        self.hbond_getter = sm.ChemTable.get_names_hbd
        self._kernel_std = _kernel_hbd
        self._kernel_alt = _kernel_hbd_fixed
        self.process_kernel()

        self.hbond_getter = sm.ChemTable.get_names_hbd_fixed
        self._kernel_std = _kernel_hbd_fixed
        self._kernel_alt = _kernel_hbd_fixed
        self.process_kernel()


    # --------------------------------------------------------------------------
    def _iter_triplets(self):
        if sm.USE_STRUCTURE_HYDROGENS:
            self.ms.system.guess_TopologyAttrs(to_guess = ["bonds"])

        for triplet in super()._iter_triplets():
            if triplet.interactor in self.processed_interactors: continue

            if sm.USE_STRUCTURE_HYDROGENS:
                for hydrogen in triplet.get_interactor_bonded_hydrogens(self.res_atoms):
                    triplet.pos_tail = triplet.pos_interactor
                    triplet.pos_head = hydrogen.position
                    self.kernel = self._kernel_alt
                    self.processed_interactors.add(triplet.interactor)
                    yield triplet

            if triplet.pos_head is None: # sm.USE_STRUCTURE_HYDROGENS falls back to "no-hydrogen" model if no hydrogens found
                self.kernel = self._kernel_std
                yield triplet


# //////////////////////////////////////////////////////////////////////////////
