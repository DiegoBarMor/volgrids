import volgrids as vg

from .stacking import SmifStacking

# //////////////////////////////////////////////////////////////////////////////
class SmifStackingPP(SmifStacking):
    """
    Stacking Probe-Probe (PP) field.

    Following overlap-gio2 approach: uses SphereKernel with radius 2.0 Å
    """

    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # PP fields use simple sphere kernel (radius 2.0A), maybe we need to recheck that or implement it on each atom instead of the whole ring
        self.kernel = vg.KernelSphere(
            radius = 2.0, # [WIP] hardcoded value -> config
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self):
        for atoms_plane in self.iter_particles():
            if len(atoms_plane) < 3: continue

            cog = atoms_plane.center_of_geometry()
            # Stamp sphere at ring center (radius 2.0 Å)
            self.kernel.stamp(self.grid, cog)


# //////////////////////////////////////////////////////////////////////////////
