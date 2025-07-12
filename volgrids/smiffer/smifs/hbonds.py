import numpy as np
from abc import ABC, abstractmethod

import volgrids.vgrids as vg
import volgrids.smiffer as sm

def _str_prev_residue(res):
    return f"segid {res.segid} and resid {res.resid - 1}"

def _str_this_residue(res):
    return f"segid {res.segid} and resid {res.resid} and resname {res.resname}"

def _str_next_residue(res):
    return f"segid {res.segid} and resid {res.resid + 1}"

def _safe_return_coords(atomgroup, sel_string):
    atoms = atomgroup.select_atoms(sel_string)
    if len(atoms) == 0: return None
    return atoms.center_of_geometry()


# //////////////////////////////////////////////////////////////////////////////
class SmifHBonds(sm.Smif, ABC):
    def populate_grid(self):
        self.kernel: vg.Kernel
        self.kernel_args: dict
        self.hbond_getter: callable
        self.init_kernel()

        self.kernel.link_to_grid(self.grid, self.ms.minCoords)
        for pos_antecedent, pos_atom_hbond in self.iter_particles():
            direction = vg.Math.normalize(pos_atom_hbond - pos_antecedent)
            self.kernel.recalculate_kernel(direction, **self.kernel_args)
            self.kernel.stamp(pos_atom_hbond, multiplication_factor = sm.ENERGY_SCALE)


    def iter_particles(self):
        atoms = self.ms.relevant_atoms
        for res in atoms.residues:
            res_atoms = atoms.select_atoms(_str_this_residue(res))
            hbond_triplets = self.hbond_getter(self.ms.chemtable, res.resname)
            if hbond_triplets is None: continue # skip weird residues

            for hbond_triplet in hbond_triplets:
                if not hbond_triplet: continue  # skip residues without HBond pairs

                name_hbond_atom = hbond_triplet[2]

                antecedent_pos = self.select_antecedent(res, res_atoms, hbond_triplet)
                sel_hbond_atom = res_atoms.select_atoms(f"name {name_hbond_atom}")
                if antecedent_pos is None: continue # skip special cases

                if not sel_hbond_atom or len(antecedent_pos) == 0: continue # skip residue artifacts
                hbond_atom_pos = sel_hbond_atom[0].position

                yield antecedent_pos, hbond_atom_pos


    @abstractmethod
    def select_antecedent(self, res, res_atoms, hbond_triplet) -> None|np.ndarray:
        return


    @abstractmethod
    def init_kernel(self):
        return


# //////////////////////////////////////////////////////////////////////////////
class SmifHBAccepts(SmifHBonds, ABC):
    def select_antecedent(self, res, res_atoms, hbond_triplet) -> None|np.ndarray:
        name_antecedent_0, name_antecedent_1, name_hbond_atom, _ = hbond_triplet
        atoms = self.ms.relevant_atoms

        ##### standard antecedents
        if not name_antecedent_1:
            return _safe_return_coords(res_atoms, f"name {name_antecedent_0}")

        ##### pseudo-antecedents
        ### special case for RNA, needs to check next residue
        if name_hbond_atom == "O3'":
            return _safe_return_coords(atoms,
                f"({_str_this_residue(res)} and name {name_antecedent_0}) or" +\
                f"({_str_next_residue(res)} and name {name_antecedent_1})"
            )

        ### other pseudo-antecedent cases
        return _safe_return_coords(res_atoms, f"name {name_antecedent_0} {name_antecedent_1}")


# //////////////////////////////////////////////////////////////////////////////
class SmifHBDonors(SmifHBonds, ABC):
    def select_antecedent(self, res, res_atoms, hbond_triplet) -> None|np.ndarray:
        name_antecedent_0, name_antecedent_1, name_hbond_atom, do_infer_H = hbond_triplet
        atoms = self.ms.relevant_atoms

        ##### pseudo-antecedents
        if name_antecedent_1:
            ### special case for terminal aminoacids
            if hbond_triplet == ("C", "CA", "N", False):
                return _safe_return_coords(atoms,
                    f"({_str_prev_residue(res)} and name C ) or " +\
                    f"({_str_this_residue(res)} and name CA)"
                )

            if do_infer_H:
                atom_ref        = res_atoms.select_atoms(f"name {name_antecedent_0}")
                atom_antecedent = res_atoms.select_atoms(f"name {name_antecedent_1}")
                atom_hbond      = res_atoms.select_atoms(f"name {name_hbond_atom  }")
                if len(atom_ref) != 1 or len(atom_antecedent) != 1 or len(atom_hbond) != 1:
                    return

                direction = atom_antecedent[0].position - atom_ref[0].position
                return atom_hbond[0].position - direction

            return _safe_return_coords(res_atoms, f"name {name_antecedent_0} {name_antecedent_1}")

        ##### standard antecedents
        ## special case for RNA, it's a donor only if there is no next residue
        if name_hbond_atom == "O3'":
            sel_next_res = atoms.select_atoms(_str_next_residue(res))
            if len(sel_next_res) > 0: return

        ## special case for RNA, it's a donor only if there is no previous residue
        if name_hbond_atom == "O5'":
            sel_prev_res = atoms.select_atoms(_str_prev_residue(res))
            if len(sel_prev_res) > 0: return

        ## other standard antecedent cases
        return _safe_return_coords(res_atoms, f"name {name_antecedent_0}")


# //////////////////////////////////////////////////////////////////////////////
class SmifHBARing(SmifHBAccepts):
    def init_kernel(self):
        kernel_radius = sm.MU_HBA[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBA
        self.kernel = vg.KernelGaussianMultivariate(kernel_radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBA, cov_inv = sm.COV_INV_HBA, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hba


# //////////////////////////////////////////////////////////////////////////////
class SmifHBDRing(SmifHBDonors):
    def init_kernel(self):
        kernel_radius = sm.MU_HBD[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD
        self.kernel = vg.KernelGaussianMultivariate(kernel_radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBD, cov_inv = sm.COV_INV_HBD, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd


# //////////////////////////////////////////////////////////////////////////////
class SmifHBDCone(SmifHBDonors):
    def init_kernel(self):
        kernel_radius = sm.MU_HBD_FIXED[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD_FIXED
        self.kernel = vg.KernelGaussianMultivariate(kernel_radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBD_FIXED, cov_inv = sm.COV_INV_HBD_FIXED, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd_fixed


# //////////////////////////////////////////////////////////////////////////////
