import volgrids as vg
import volgrids.smiffer as sm

from ._core.hydro import SmifHydro

# //////////////////////////////////////////////////////////////////////////////
class SmifHydrophilic(SmifHydro):
    def populate_grid(self, grid: vg.Grid) -> None:
        grid.reset()
        radius = sm.PARAM_HPHIL_DIST_MU + sm.MISC_KERNEL_GAUSSIAN_SIGMAS * sm.PARAM_HPHIL_DIST_SIGMA
        kernel = vg.KernelGaussianUnivariateDist(
            radius, self.ms.get_deltas(), vg.FLOAT_DTYPE, sm.PARAMS_HPHIL
        )
        for particle, mul_factor in self.iter_particles():
            if mul_factor > 0: continue
            kernel.stamp(grid, particle.position, multiply_by = -mul_factor)


# //////////////////////////////////////////////////////////////////////////////
