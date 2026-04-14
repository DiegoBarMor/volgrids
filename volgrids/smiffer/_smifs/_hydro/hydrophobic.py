import volgrids as vg
import volgrids.smiffer as sm

from .hydro import SmifHydro

# //////////////////////////////////////////////////////////////////////////////
class SmifHydrophobic(SmifHydro):
    def populate_grid(self):
        radius = sm.MU_HYDROPHOBIC + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_HYDROPHOBIC
        kernel = vg.KernelGaussianUnivariateDist(
            radius, self.ms.get_deltas(), vg.FLOAT_DTYPE, sm.PARAMS_HPHOB
        )
        for particle, mul_factor in self.iter_particles():
            if mul_factor < 0: continue
            kernel.stamp(self.grid, particle.position, multiply_by =  mul_factor)


# //////////////////////////////////////////////////////////////////////////////
