import volpot
import numpy as np

from settings import HPHOB_RNA_PHOSPHATE, HPHOB_RNA_SUGAR, \
    MU_HYDROPHILIC, SIGMA_HYDROPHILIC, MU_HYDROPHOBIC, SIGMA_HYDROPHOBIC, GAUSSIAN_KERNEL_SIGMAS


# //////////////////////////////////////////////////////////////////////////////
class SPG_Hydro(volpot.StatisticalPotentialGrid):
    def iter_particles(self):
        for atom in self.ms.relevant_atoms:
            # atom is part of a protein
            if not self.ms.isNucleic:
                mul_factor = volpot.ww_scale[atom.resname]

            # atom is in the nitrogenous base (RNA)
            elif atom.name in volpot.nucleic_bases[atom.resname]:
                mul_factor = volpot.ww_scale[atom.resname]

            # atom is part of the phosphate backbone (RNA)
            elif atom.name in volpot.nucleic_backbone_phosphate:
                mul_factor = HPHOB_RNA_PHOSPHATE

            # atom is part of the sugar (RNA)
            elif atom.name in volpot.nucleic_backbone_sugar:
                mul_factor = HPHOB_RNA_SUGAR

            else: continue # unknown atom type

            yield atom, mul_factor


# //////////////////////////////////////////////////////////////////////////////
class SPG_Hydrophilic(SPG_Hydro):
    POTENTIAL_TYPE = "hydrophilic"

    def populate_grid(self):
        radius_hphil = MU_HYDROPHILIC + GAUSSIAN_KERNEL_SIGMAS * SIGMA_HYDROPHILIC
        gk_hphil = volpot.GaussianKernel(MU_HYDROPHILIC, SIGMA_HYDROPHILIC, radius_hphil, self.ms.deltas, np.float32)
        gk_hphil.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor > 0: continue
            gk_hphil.stamp(particle.position, multiplication_factor = -mul_factor)


# //////////////////////////////////////////////////////////////////////////////
class SPG_Hydrophobic(SPG_Hydro):
    POTENTIAL_TYPE = "hydrophobic"

    def populate_grid(self):
        radius_hphob = MU_HYDROPHOBIC + GAUSSIAN_KERNEL_SIGMAS * SIGMA_HYDROPHOBIC
        gk_hphob = volpot.GaussianKernel(MU_HYDROPHOBIC, SIGMA_HYDROPHOBIC, radius_hphob, self.ms.deltas, np.float32)
        gk_hphob.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor < 0: continue
            gk_hphob.stamp(particle.position, multiplication_factor = mul_factor)


# //////////////////////////////////////////////////////////////////////////////
