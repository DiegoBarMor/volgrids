import volgrids as vg
import volgrids.smiffer as sm

from .hbonds import SmifHBonds
from .triplet import Triplet

# //////////////////////////////////////////////////////////////////////////////
class SmifHBDonorsPP(SmifHBonds):
    """
    Hydrogen Bond Donor Probe-Probe (PP) field.

    Following overlap-gio2 approach: uses SphereKernel with radius 2.0 Å, to be rechecked
    """

    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # PP fields use simple sphere kernel (radius 2.0 Å)
        self.kernel = vg.KernelSphere(
            radius = 2.0, # [WIP] hardcoded value -> config
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )
        self.hbond_getter = sm.ParserChemTable.get_names_hbd

    # --------------------------------------------------------------------------
    def find_tail_head_positions(self, triplet: Triplet) -> None:
        # PP fields only use head position (donor site)
        triplet.set_pos_head(self.res_atoms)
        # No tail position needed for spherical PP fieldsx

    # --------------------------------------------------------------------------
    def populate_grid(self):
        """Populate grid with spherical accessibility regions."""
        for triplet in self._iter_triplets():
            if triplet.pos_interactor is None: continue
            # Stamp sphere at donor position (radius 2.0 Å)
            self.kernel.stamp(self.grid, triplet.pos_interactor)



# //////////////////////////////////////////////////////////////////////////////
