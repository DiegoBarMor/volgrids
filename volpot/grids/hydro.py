import numpy as np
import volpot as vp

# //////////////////////////////////////////////////////////////////////////////
class GridHydro(vp.GridSMIF):
    def iter_particles(self):
        for atom in self.ms.get_relevant_atoms():
            # atom is part of a protein
            if not self.ms.isNucleic:
                mul_factor = vp.ww_scale[atom.resname]

            # atom is in the nitrogenous base (RNA)
            elif atom.name in vp.nucleic_bases[atom.resname]:
                mul_factor = vp.ww_scale[atom.resname]

            # atom is part of the phosphate backbone (RNA)
            elif atom.name in vp.nucleic_backbone_phosphate:
                mul_factor = vp.HPHOB_RNA_PHOSPHATE

            # atom is part of the sugar (RNA)
            elif atom.name in vp.nucleic_backbone_sugar:
                mul_factor = vp.HPHOB_RNA_SUGAR

            else: continue # unknown atom type

            yield atom, mul_factor


# //////////////////////////////////////////////////////////////////////////////
class GridHydrophilic(GridHydro):
    def get_type(self):
        return "hydrophilic"

    def populate_grid(self):
        radius_hphil = vp.MU_HYDROPHILIC + vp.GAUSSIAN_KERNEL_SIGMAS * vp.SIGMA_HYDROPHILIC
        gk_hphil = vp.KernelGaussianUnivariate(
            vp.MU_HYDROPHILIC, vp.SIGMA_HYDROPHILIC,
            radius_hphil, self.ms.deltas, vp.FLOAT_DTYPE
        )
        gk_hphil.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor > 0: continue
            gk_hphil.stamp(particle.position, multiplication_factor = -mul_factor)


# //////////////////////////////////////////////////////////////////////////////
class GridHydrophobic(GridHydro):
    def get_type(self):
        return "hydrophobic"

    def populate_grid(self):
        radius_hphob = vp.MU_HYDROPHOBIC + vp.GAUSSIAN_KERNEL_SIGMAS * vp.SIGMA_HYDROPHOBIC
        gk_hphob = vp.KernelGaussianUnivariate(
            vp.MU_HYDROPHOBIC, vp.SIGMA_HYDROPHOBIC,
            radius_hphob, self.ms.deltas, vp.FLOAT_DTYPE
        )
        gk_hphob.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor < 0: continue
            gk_hphob.stamp(particle.position, multiplication_factor = mul_factor)


# //////////////////////////////////////////////////////////////////////////////
