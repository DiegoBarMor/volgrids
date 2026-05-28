import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class OgHydrophobic(sm.SmifHydrophobic):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel = vg.KernelSphere(
            radius = su.OG_RADIUS_HYDROPHOBIC,
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid):
        """Populate grid with spherical accessibility regions."""
        grid.reset()
        for atom, logp_value in self.iter_particles():
            ### [TODO] check this
            if logp_value <= 0: continue  # Only hydrophobic atoms
            self.kernel.stamp(grid, atom.position)


# //////////////////////////////////////////////////////////////////////////////
