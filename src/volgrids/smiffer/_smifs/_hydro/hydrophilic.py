import volgrids as vg
import volgrids.smiffer as sm

from .hydro import SmifHydro

# //////////////////////////////////////////////////////////////////////////////
class SmifHydrophilic(SmifHydro):
    def populate_grid(self):
        radius = sm.MU_HYDROPHILIC + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_HYDROPHILIC
        kernel = vg.KernelGaussianUnivariateDist(
            radius, self.ms.deltas, vg.FLOAT_DTYPE, sm.PARAMS_HPHIL
        )
        kernel.link_to_grid(self.grid, self.ms.minCoords)

        for particle, mul_factor in self.iter_particles():
            if mul_factor > 0: continue
            kernel.stamp(particle.position, multiplication_factor = -mul_factor)


# //////////////////////////////////////////////////////////////////////////////
