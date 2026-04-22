import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class OgHydrophobic(sm.SmifHydrophobic):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel = vg.KernelSphere(
            radius = 2.0, # [TODO] hardcoded value -> config
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self):
        """Populate grid with spherical accessibility regions."""
        for atom, logp_value in self.iter_particles():
            ### [TODO] check this
            if logp_value <= 0: continue  # Only hydrophobic atoms
            self.kernel.stamp(self.grid, atom.position)


# //////////////////////////////////////////////////////////////////////////////
