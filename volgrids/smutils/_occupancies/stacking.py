import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class OgStacking(sm.SmifStacking):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kernel = vg.KernelSphere(
            radius = su.OG_RADIUS_STACKING,
            deltas = self.ms.get_deltas(),
            dtype = vg.FLOAT_DTYPE
        )

    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid) -> None:
        grid.reset()
        for atoms_plane in self.iter_particles():
            for atom in atoms_plane:
                self.kernel.stamp(grid, atom.position)


# //////////////////////////////////////////////////////////////////////////////
