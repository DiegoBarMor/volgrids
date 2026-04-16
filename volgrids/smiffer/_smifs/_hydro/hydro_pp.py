import volgrids as vg
import volgrids.smiffer as sm

from .hydrophobic import SmifHydrophobic

# //////////////////////////////////////////////////////////////////////////////
class SmifHydroPP(SmifHydrophobic):
    """
    Hydrophobic Probe-Probe (PP) field.

    Following overlap-gio2 approach: uses SphereKernel with radius 2.0 Å
    """

    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # PP fields use simple sphere kernel (radius 2.0 Å)
        self.kernel = vg.KernelSphere(
            radius = 2.0,
            deltas = self.ms.deltas,
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self):
        """Populate grid with spherical accessibility regions."""
        for atom, logp_value in self.iter_particles():
            # NOT SURE IT'S GOOD, PRETTY SURE IT'S BAD
            if logp_value > 0:  # Only hydrophobic atoms
                # Stamp sphere at atom position (radius 2.0 Å)
                self.kernel.stamp(self, atom.position)

    # --------------------------------------------------------------------------
    def save_data(self, folder_out, title):
        super().save_data(folder_out, title)


# //////////////////////////////////////////////////////////////////////////////
