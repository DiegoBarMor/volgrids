import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class OgHBAccepts(sm._smifs_core.SmifHBonds):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel = vg.KernelSphere(
            radius = su.RADIUS_OCCUPANCY_OG,
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )
        self.hbond_getter = sm.ParserChemTable.get_names_hba

    # --------------------------------------------------------------------------
    def find_tail_head_positions(self, triplet: sm._smifs_core.Triplet) -> None:
        ### OGs only use head position (acceptor site) --> no tail position needed
        triplet.set_pos_head(self.res_atoms)

    # --------------------------------------------------------------------------
    def populate_grid(self):
        """Populate grid with spherical accessibility regions."""
        for triplet in self._iter_triplets():
            if triplet.pos_interactor is None: continue
            self.kernel.stamp(self.grid, triplet.pos_interactor)


# //////////////////////////////////////////////////////////////////////////////
