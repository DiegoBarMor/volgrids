import numpy as np
import volgrids as vg
import volgrids.smiffer as sm
from abc import ABC, abstractmethod

# //////////////////////////////////////////////////////////////////////////////
class GridHBonds(vg.Grid, ABC):
    def populate_grid(self):
        self.kernel: vg.Kernel
        self.kernel_args: dict
        self.hbond_getter: callable
        self.init_kernel()

        self.kernel.link_to_grid(self.grid, self.ms.minCoords)
        for pos_antecedent, pos_atom_hbond in self.iter_particles():
            direction = vg.Math.normalize(pos_atom_hbond - pos_antecedent)
            self.kernel.recalculate_kernel(direction, **self.kernel_args)
            self.kernel.stamp(pos_atom_hbond)


    def iter_particles(self):
        atoms = self.ms.relevant_atoms
        for res in atoms.residues:
            res_atoms = atoms.select_atoms(f"segid {res.segid} and resid {res.resid} and resname {res.resname}")
            hbond_triplets = self.hbond_getter(self.ms.chemtable, res.resname)
            if hbond_triplets is None: continue # skip weird residues

            for hbond_triplet in hbond_triplets:
                if not hbond_triplet: continue  # skip residues without HBond pairs

                name_hbond_atom = hbond_triplet[2]

                sel_antecedent = self.select_antecedent(res, res_atoms, hbond_triplet)
                sel_hbond_atom = res_atoms.select_atoms(f"name {name_hbond_atom}")
                if sel_antecedent is None: continue # skip special cases

                if not sel_hbond_atom or not sel_antecedent: continue # skip residue artifacts
                antecedent_pos = np.mean(sel_antecedent.positions, axis = 0)
                hbond_atom_pos = sel_hbond_atom[0].position

                yield antecedent_pos, hbond_atom_pos


    @abstractmethod
    def select_antecedent(self, res, res_atoms, hbond_triplet):
        return


    @abstractmethod
    def init_kernel(self):
        return


# //////////////////////////////////////////////////////////////////////////////
class GridHBAccepts(GridHBonds, ABC):
    def select_antecedent(self, res, res_atoms, hbond_triplet):
        name_antecedent_0, name_antecedent_1, name_hbond_atom = hbond_triplet
        atoms = self.ms.relevant_atoms

        ##### standard antecedents
        if not name_antecedent_1:
            return res_atoms.select_atoms(f"name {name_antecedent_0}")

        ##### pseudo-antecedents
        ### special case for RNA, needs to check next residue
        if name_hbond_atom == "O3'":
            return atoms.select_atoms(
                f"(segid {res.segid} and resid {res.resid} and resname {res.resname} and name {name_antecedent_0}) or" +\
                f"(segid {res.segid} and resid {res.resid + 1} " +                 f"and name {name_antecedent_1})"
            )

        ### other pseudo-antecedent cases
        return res_atoms.select_atoms(f"name {name_antecedent_0} {name_antecedent_1}")


# //////////////////////////////////////////////////////////////////////////////
class GridHBDonors(GridHBonds, ABC):
    def select_antecedent(self, res, res_atoms, hbond_triplet):
        name_antecedent_0, name_antecedent_1, name_hbond_atom = hbond_triplet
        atoms = self.ms.relevant_atoms

        ##### pseudo-antecedents
        if name_antecedent_1:
            return res_atoms.select_atoms(f"name {name_antecedent_0} {name_antecedent_1}")

        ##### standard antecedents
        ## special case for RNA, it's a donor only if there is no next residue
        if name_hbond_atom == "O3'":
            sel_next_res = atoms.select_atoms(
                f"segid {res.segid} and resid {res.resid + 1}"
            )
            if len(sel_next_res) > 0: return

        ## special case for RNA, it's a donor only if there is no previous residue
        if name_hbond_atom == "O5'":
            sel_prev_res = atoms.select_atoms(
                f"segid {res.segid} and resid {res.resid - 1}"
            )
            if len(sel_prev_res) > 0: return

        ## other standard antecedent cases
        return res_atoms.select_atoms(f"name {name_antecedent_0}")


# //////////////////////////////////////////////////////////////////////////////
class GridHBARing(GridHBAccepts):
    def init_kernel(self):
        kernel_radius = sm.MU_HBA[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBA
        self.kernel = vg.KernelGaussianMultivariate(kernel_radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBA, cov_inv = sm.COV_INV_HBA, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hba


# //////////////////////////////////////////////////////////////////////////////
class GridHBDRing(GridHBDonors):
    def init_kernel(self):
        kernel_radius = sm.MU_HBD[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD
        self.kernel = vg.KernelGaussianMultivariate(kernel_radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBD, cov_inv = sm.COV_INV_HBD, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd


# //////////////////////////////////////////////////////////////////////////////
class GridHBDCone(GridHBDonors):
    def init_kernel(self):
        kernel_radius = sm.MU_HBD_FIXED[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBD_FIXED
        self.kernel = vg.KernelGaussianMultivariate(kernel_radius, self.ms.deltas, vg.FLOAT_DTYPE)

        self.kernel_args = dict(mu = sm.MU_HBD_FIXED, cov_inv = sm.COV_INV_HBD_FIXED, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd_fixed


# //////////////////////////////////////////////////////////////////////////////
