import numpy as np
import volgrids.vgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class GridHBonds(vg.Grid):
    def populate_grid(self):
        self.radius: float
        self.hbond_getter: callable
        self.kernel_args: dict

        gk = vg.KernelGaussianMultivariate(self.radius, self.ms.deltas, vg.FLOAT_DTYPE)
        gk.link_to_grid(self.grid, self.ms.minCoords)
        for pos_antecedent, pos_atom_hbond in self.iter_particles():
            direction = vg.normalize(pos_atom_hbond - pos_antecedent)
            gk.recalculate_kernel(direction, **self.kernel_args)
            gk.stamp(pos_atom_hbond)


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


    def select_antecedent(self, res, res_atoms, hbond_triplet):
        return # override to use


# //////////////////////////////////////////////////////////////////////////////
class GridHBAccepts(GridHBonds):
    def get_type(self):
        return "hbacceptors"

    def populate_grid(self):
        self.radius = sm.MU_HBA[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBOND
        self.kernel_args = dict(mu = sm.MU_HBA, cov_inv = sm.COV_INV_HBA, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hba
        super().populate_grid()

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
class GridHBDonors(GridHBonds):
    def get_type(self):
        return "hbdonors"

    def populate_grid(self):
        self.radius = sm.MU_HBD[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_HBOND
        self.kernel_args = dict(mu = sm.MU_HBD, cov_inv = sm.COV_INV_HBD, isStacking = False)
        self.hbond_getter = sm.ChemTable.get_names_hbd
        super().populate_grid()

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
