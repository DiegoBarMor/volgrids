import numpy as np
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class GridHydro(vg.GridSMIF):
    def iter_particles(self):
        for atom in self.ms.get_relevant_atoms():
            # atom is part of a protein
            if not self.ms.isNucleic:
                mul_factor = vg.ww_scale[atom.resname]

            # atom is in the nitrogenous base (RNA)
            elif atom.name in vg.nucleic_bases[atom.resname]:
                mul_factor = vg.ww_scale[atom.resname]

            # atom is part of the phosphate backbone (RNA)
            elif atom.name in vg.nucleic_backbone_phosphate:
                mul_factor = vg.HPHOB_RNA_PHOSPHATE

            # atom is part of the sugar (RNA)
            elif atom.name in vg.nucleic_backbone_sugar:
                mul_factor = vg.HPHOB_RNA_SUGAR

            else: continue # unknown atom type

            yield atom, mul_factor


# //////////////////////////////////////////////////////////////////////////////
class GridHydrophilic(GridHydro):
    def get_type(self):
        return "hydrophilic"

    def populate_grid(self):
        radius_hphil = vg.MU_HYDROPHILIC + vg.GAUSSIAN_KERNEL_SIGMAS * vg.SIGMA_HYDROPHILIC
        gk_hphil = vg.KernelGaussianUnivariate(
            vg.MU_HYDROPHILIC, vg.SIGMA_HYDROPHILIC,
            radius_hphil, self.ms.deltas, vg.FLOAT_DTYPE
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
        radius_hphob = vg.MU_HYDROPHOBIC + vg.GAUSSIAN_KERNEL_SIGMAS * vg.SIGMA_HYDROPHOBIC
        gk_hphob = vg.KernelGaussianUnivariate(
            vg.MU_HYDROPHOBIC, vg.SIGMA_HYDROPHOBIC,
            radius_hphob, self.ms.deltas, vg.FLOAT_DTYPE
        )
        gk_hphob.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor < 0: continue
            gk_hphob.stamp(particle.position, multiplication_factor = mul_factor)


# //////////////////////////////////////////////////////////////////////////////
