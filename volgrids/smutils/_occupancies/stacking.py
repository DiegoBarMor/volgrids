import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class OgStacking(sm.SmifStacking):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        ### [TODO] sphere on each atom instead of the whole ring
        self.kernel = vg.KernelSphere(
            radius = 2.0, # [TODO] hardcoded value -> config
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self):
        for atoms_plane in self.iter_particles():
            if len(atoms_plane) < 3: continue

            cog = atoms_plane.center_of_geometry()
            ### [WIP] Stamp sphere at ring center (radius 2.0 Å)
            self.kernel.stamp(self.grid, cog)


# //////////////////////////////////////////////////////////////////////////////
