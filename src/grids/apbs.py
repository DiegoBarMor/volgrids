import numpy as np
from gridData import Grid

from src.utilities import Timer, interpolate_3d
from src.grids.potential_grids import PotentialGrid

from settings import APBS_MIN_CUTOFF, APBS_MAX_CUTOFF, DO_LOG_APBS


# //////////////////////////////////////////////////////////////////////////////
class PPG_APBS(PotentialGrid):
    POTENTIAL_TYPE = "apbs"

    def __init__(self, ms):
        if str(ms.path_apbs) == '.':
            return print("...>>> APBS output not provided. Electrostatic potential skipped.", flush = True)

        if not ms.path_apbs.exists():
            return print(f"...XXX APBS output not found, '{ms.path_apbs}' doesn't exist. Electrostatic potential skipped.", flush = True)

        super().__init__(ms)

        if DO_LOG_APBS:
            timer = Timer("...>>> Creating APBS-LOG potential grid...", flush = True)
            self.apply_logabs_transform()
            timer.end()
        else:
            print("...>>> APBS-LOG potential grid skipped.", flush = True)


    def populate_grid(self):
        apbs = Grid(self.ms.path_apbs)

        xmin0, ymin0, zmin0 = apbs.origin
        xmax0, ymax0, zmax0 = apbs.origin + apbs.delta * apbs.grid.shape
        xres0, yres0, zres0 = apbs.grid.shape

        xmin1, ymin1, zmin1 = self.ms.minCoords
        xmax1, ymax1, zmax1 = self.ms.maxCoords
        xres1, yres1, zres1 = self.ms.resolution

        self.grid = interpolate_3d(
            x0 = np.linspace(xmin0, xmax0, xres0),
            y0 = np.linspace(ymin0, ymax0, yres0),
            z0 = np.linspace(zmin0, zmax0, zres0),
            data_0 = apbs.grid,
            new_coords = np.mgrid[
                xmin1 : xmax1 : complex(0, xres1),
                ymin1 : ymax1 : complex(0, yres1),
                zmin1 : zmax1 : complex(0, zres1),
            ].T
        ).astype(np.float32)


    def apply_logabs_transform(self):
        logpos = np.log10( self.grid[self.grid > 0])
        logneg = np.log10(-self.grid[self.grid < 0])

        ##### APPLY CUTOFFS
        logpos[logpos < APBS_MIN_CUTOFF] = APBS_MIN_CUTOFF
        logneg[logneg < APBS_MIN_CUTOFF] = APBS_MIN_CUTOFF
        logpos[logpos > APBS_MAX_CUTOFF] = APBS_MAX_CUTOFF
        logneg[logneg > APBS_MAX_CUTOFF] = APBS_MAX_CUTOFF

        ##### SHIFT VALUES TO 0
        logpos -= APBS_MIN_CUTOFF
        logneg -= APBS_MIN_CUTOFF

        ##### REVERSE SIGN OF LOG(ABS(GRID_NEG)) AND DOUBLE BOTH
        logpos *=  2 # this way the range of points varies between
        logneg *= -2 # 2*APBS_MIN_CUTOFF and 2*APBS_MAX_CUTOFF

        ##### RESULT
        self.grid[self.grid > 0] = logpos
        self.grid[self.grid < 0] = logneg

        ##### additional adjustments before saving to json
        self.pack_data()
        self.data["minPotential"] = float(-2 * (APBS_MAX_CUTOFF - APBS_MIN_CUTOFF))
        self.data["maxPotential"] = float( 2 * (APBS_MAX_CUTOFF - APBS_MIN_CUTOFF))
        self.save_data(override_prefix = "apbslog")


# //////////////////////////////////////////////////////////////////////////////
