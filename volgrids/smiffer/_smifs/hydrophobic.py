import volgrids as vg
import volgrids.smiffer as smf

from ._core.hydro import SmifHydro

# //////////////////////////////////////////////////////////////////////////////
class SmifHydrophobic(SmifHydro):
    def populate_grid(self, grid: vg.Grid) -> None:
        grid.reset()
        radius = vg.CFG.param_hphob_dist_mu + vg.CFG.misc_kernel_gaussian_sigmas * vg.CFG.param_hbhob_dist_sigma
        kernel = vg.KernelGaussianUnivariateDist(
            radius, self.ms.get_deltas(), vg.FLOAT_DTYPE, smf.PARAMS_HPHOB
        )
        for particle, mul_factor in self.iter_particles():
            if mul_factor < 0: continue
            kernel.stamp(grid, particle.position, multiply_by =  mul_factor)


# //////////////////////////////////////////////////////////////////////////////
