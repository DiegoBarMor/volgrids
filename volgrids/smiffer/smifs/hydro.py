from abc import ABC

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifHydro(sm.Smif, ABC):
    def iter_particles(self):
        for atom in self.ms.relevant_atoms:
            factor_res  = self.ms.chemtable.get_residue_hphob(atom)
            factor_atom = self.ms.chemtable.get_atom_hphob(atom)

            if (factor_res is None) and (factor_atom is None):
                continue # skip atoms with unknown name and resname

            if factor_res  is None: factor_res  = 1
            if factor_atom is None: factor_atom = 1

            yield atom, factor_res * factor_atom #/ len(atom.residue.atoms)


# //////////////////////////////////////////////////////////////////////////////
class SmifHydrophilic(SmifHydro):
    def populate_grid(self):
        radius_hphil = sm.MU_HYDROPHILIC + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_HYDROPHILIC
        gk_hphil = vg.KernelGaussianUnivariate(
            sm.MU_HYDROPHILIC, sm.SIGMA_HYDROPHILIC,
            radius_hphil, self.ms.deltas, vg.FLOAT_DTYPE
        )
        gk_hphil.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor > 0: continue
            gk_hphil.stamp(particle.position, multiplication_factor = -mul_factor)


# //////////////////////////////////////////////////////////////////////////////
class SmifHydrophobic(SmifHydro):
    def populate_grid(self):
        radius_hphob = sm.MU_HYDROPHOBIC + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_HYDROPHOBIC
        gk_hphob = vg.KernelGaussianUnivariate(
            sm.MU_HYDROPHOBIC, sm.SIGMA_HYDROPHOBIC,
            radius_hphob, self.ms.deltas, vg.FLOAT_DTYPE
        )
        gk_hphob.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor < 0: continue
            gk_hphob.stamp(particle.position, multiplication_factor = mul_factor)


# //////////////////////////////////////////////////////////////////////////////
