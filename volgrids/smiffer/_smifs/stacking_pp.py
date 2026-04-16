import volgrids as vg
import volgrids.smiffer as sm

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
            radius = 2.0,
            deltas = self.ms.deltas,
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self):
        for ring_atoms in self.iter_particles():
            if len(ring_atoms) >= 3:
                # Get center of geometry of aromatic ring
                center_of_geometry = ring_atoms.center_of_geometry()

                # Stamp sphere at ring center (radius 2.0 Å)
                self.kernel.stamp(self, center_of_geometry)

    # --------------------------------------------------------------------------
    def save_data(self, folder_out, title):
        super().save_data(folder_out, title)


# //////////////////////////////////////////////////////////////////////////////
