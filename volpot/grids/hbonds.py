import numpy as np
import volpot as vp

# //////////////////////////////////////////////////////////////////////////////
class SPG_HBonds(vp.StatisticalPotentialGrid):
    def populate_grid(self):
        self.radius: float
        self.table_hbond: dict
        self.kernel_args: dict

        gk = vp.MultiGaussianKernel(self.radius, self.ms.deltas, np.float32)
        gk.link_to_grid(self.grid, self.ms.minCoords)
        for pos_antecedent, pos_atom_hbond in self.iter_particles():
            direction = vp.normalize(pos_atom_hbond - pos_antecedent)
            gk.recalculate_kernel(direction, **self.kernel_args)
            gk.stamp(pos_atom_hbond)


    def iter_particles(self):
        for res in self.ms.relevant_atoms.residues:
            res_atoms = self.ms.relevant_atoms.select_atoms(f"segid {res.segid} and resid {res.resid} and resname {res.resname}")
            hbond_pairs = self.table_hbond.get(res.resname)
            if hbond_pairs is None: continue # skip weird residues

            for hbond_pair in hbond_pairs:
                if not hbond_pair: continue  # skip residues without HBond pairs
                name_antecedent, name_hbond_atom = hbond_pair

                sel_antecedent = self.select_antecedent(res, res_atoms, name_antecedent, name_hbond_atom)
                sel_hbond_atom = res_atoms.select_atoms(f"name {name_hbond_atom}")
                if sel_antecedent is None: continue # skip special cases

                if not sel_hbond_atom or not sel_antecedent: continue # skip residue artifacts
                antecedent_pos = np.mean(sel_antecedent.positions, axis = 0)
                hbond_atom_pos = sel_hbond_atom[0].position

                yield antecedent_pos, hbond_atom_pos


    def select_antecedent(self, res, res_atoms, name_antecedent, name_hbond_atom):
        return # override to use


# //////////////////////////////////////////////////////////////////////////////
class SPG_HB_Accepts(SPG_HBonds):
    POTENTIAL_TYPE = "hbacceptors"

    def populate_grid(self):
        self.radius = vp.MU_HBA[1] + vp.GAUSSIAN_KERNEL_SIGMAS * vp.SIGMA_DIST_HBOND
        self.table_hbond = vp.rna_hba if self.ms.isNucleic else vp.prot_hba
        self.kernel_args = dict(mu = vp.MU_HBA, cov_inv = vp.COV_INV_HBA, isStacking = False)
        super().populate_grid()

    def select_antecedent(self, res, res_atoms, name_antecedent, name_hbond_atom):
        assert isinstance(name_antecedent, tuple)

        ##### standard antecedents
        if len(name_antecedent) == 1:
            return res_atoms.select_atoms(f"name {name_antecedent[0]}")

        ##### pseudo-antecedents
        if len(name_antecedent) == 2:
            name_antecedent_0, name_antecedent_1 = name_antecedent

            ### special case for RNA, needs to check next residue
            if name_hbond_atom == "O3'":
                return self.ms.relevant_atoms.select_atoms(
                    f"(segid {res.segid} and resid {res.resid} and resname {res.resname} and name {name_antecedent_0}) or" +\
                    f"(segid {res.segid} and resid {res.resid + 1} " +                 f"and name {name_antecedent_1})"
                )

            ### other pseudo-antecedent cases
            return res_atoms.select_atoms(f"name {name_antecedent_0} {name_antecedent_1}")

        raise ValueError(f"Invalid number of antecent atoms: {name_antecedent}")


# //////////////////////////////////////////////////////////////////////////////
class SPG_HB_Donors(SPG_HBonds):
    POTENTIAL_TYPE = "hbdonors"

    def populate_grid(self):
        self.radius = vp.MU_HBD[1] + vp.GAUSSIAN_KERNEL_SIGMAS * vp.SIGMA_DIST_HBOND
        self.table_hbond = vp.rna_hbd if self.ms.isNucleic else vp.prot_hbd
        self.kernel_args = dict(mu = vp.MU_HBD, cov_inv = vp.COV_INV_HBD, isStacking = False)
        super().populate_grid()

    def select_antecedent(self, res, res_atoms, name_antecedent, name_hbond_atom):
        assert isinstance(name_antecedent, tuple)

        ##### pseudo-antecedents
        if len(name_antecedent) == 2:
            name_antecedent_0, name_antecedent_1 = name_antecedent
            return res_atoms.select_atoms(f"name {name_antecedent_0} {name_antecedent_1}")

        ##### standard antecedents
        if len(name_antecedent) != 1:
            raise ValueError(f"Invalid number of antecent atoms: {name_antecedent}")

        ## special case for RNA, it's a donor only if there is no next residue
        if name_hbond_atom == "O3'":
            sel_next_res = self.ms.relevant_atoms.select_atoms(
                f"segid {res.segid} and resid {res.resid + 1}"
            )
            if len(sel_next_res) > 0: return

        ## special case for RNA, it's a donor only if there is no previous residue
        if name_hbond_atom == "O5'":
            sel_prev_res = self.ms.relevant_atoms.select_atoms(
                f"segid {res.segid} and resid {res.resid - 1}"
            )
            if len(sel_prev_res) > 0: return

        ## other standard antecedent cases
        return res_atoms.select_atoms(f"name {name_antecedent[0]}")


# //////////////////////////////////////////////////////////////////////////////
