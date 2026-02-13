import numpy as np
import pandas as pd

import volgrids as vg
import volgrids.veins as ve

# //////////////////////////////////////////////////////////////////////////////
class GridVolumetricForce(vg.Grid):
    RADIUS_FIX = 0.4 # [WIP] necessary?
    WIDTH_CYLINDER = 0.25

    # --------------------------------------------------------------------------
    def __init__(self, ms: vg.MolSystem, coords: np.ndarray, forces: np.ndarray):
        assert coords.shape == forces.shape
        super().__init__(ms)
        self.coords = coords
        self.forces = forces
        self.magnitudes = np.linalg.norm(forces, axis = 1)

        max_mag = np.max(self.magnitudes)
        if max_mag == 0.0:
            print(">>> WARNING: All forces have zero magnitude. The grid will be empty.", flush = True)

        normalization_factor = 1/max_mag if max_mag > 0 else 0.0
        self.forces *= normalization_factor * ve.FORCE_VECTOR_LEN
        self.magnitudes *= normalization_factor * ve.FORCE_VECTOR_LEN


    # --------------------------------------------------------------------------
    def populate_grid(self):
        for pos,force,mag in zip(self.coords, self.forces, self.magnitudes):
            if mag < ve.FORCE_CUTOFF: continue
            kernel = vg.KernelCylinder(
                length = mag * self.RADIUS_FIX,
                direction = force,
                width = self.WIDTH_CYLINDER,
                deltas = self.ms.deltas, dtype = np.float32
            )
            kernel.link_to_grid(self.grid, self.ms.minCoords)
            kernel.stamp(
                center_stamp_at = pos + force / 2,
                multiplication_factor = mag,
                operation = "max"
            )


# //////////////////////////////////////////////////////////////////////////////
