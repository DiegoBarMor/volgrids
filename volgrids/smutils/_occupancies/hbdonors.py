import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class OgHBDonors(sm._smifs_core.SmifHBonds):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel = vg.KernelSphere(
            radius = su.OG_RADIUS_HBD,
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )
        self.hbond_getter = sm.ParserChemTable.get_names_hbd

    # --------------------------------------------------------------------------
    def can_be_interactor(self, triplet) -> bool:
        return sm.SmifHBDonors.can_be_interactor(self, triplet)

    # --------------------------------------------------------------------------
    def find_tail_head_positions(self, triplet: sm._smifs_core.Triplet) -> None:
        ### OGs only use head position (donor site) --> no tail position needed
        triplet.set_pos_head(self.res_atoms)

    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid) -> None:
        """Populate grid with spherical accessibility regions."""
        grid.reset()
        for triplet in self._iter_triplets():
            if triplet.pos_interactor is None: continue
            self.kernel.stamp(grid, triplet.pos_interactor)


# //////////////////////////////////////////////////////////////////////////////
