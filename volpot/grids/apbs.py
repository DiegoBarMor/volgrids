import numpy as np
import volpot as vp
from gridData import Grid

# //////////////////////////////////////////////////////////////////////////////
class GridAPBS(vp.Grid):
    def get_type(self):
        return "apbs"

    def run(self, trimmer: "vp.GridTrimmer"):
        if str(self.ms.path_apbs) == '.':
            print("...--- APBS output not provided. Electrostatic potential skipped.", flush = True)
            return

        if not self.ms.path_apbs.exists():
            print(f"...XXX APBS output not found, '{self.ms.path_apbs}' doesn't exist. Electrostatic potential skipped.", flush = True)
            return

        super().run(trimmer)

    def populate_grid(self):
        apbs = Grid(self.ms.path_apbs)

        xmin0, ymin0, zmin0 = apbs.origin
        xmax0, ymax0, zmax0 = apbs.origin + apbs.delta * apbs.grid.shape
        xres0, yres0, zres0 = apbs.grid.shape

        xmin1, ymin1, zmin1 = self.ms.minCoords
        xmax1, ymax1, zmax1 = self.ms.maxCoords
        xres1, yres1, zres1 = self.ms.resolution

        self.grid = vp.interpolate_3d(
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
        if self.is_empty():
            print(f"...--- APBS potential grid is empty. Skipping logabs transform.", flush = True)
            return

        timer = vp.Timer("...>>> Creating APBS-LOG potential grid...")
        logpos = np.log10( self.grid[self.grid > 0])
        logneg = np.log10(-self.grid[self.grid < 0])

        ##### APPLY CUTOFFS
        logpos[logpos < vp.APBS_MIN_CUTOFF] = vp.APBS_MIN_CUTOFF
        logneg[logneg < vp.APBS_MIN_CUTOFF] = vp.APBS_MIN_CUTOFF
        logpos[logpos > vp.APBS_MAX_CUTOFF] = vp.APBS_MAX_CUTOFF
        logneg[logneg > vp.APBS_MAX_CUTOFF] = vp.APBS_MAX_CUTOFF

        ##### SHIFT VALUES TO 0
        logpos -= vp.APBS_MIN_CUTOFF
        logneg -= vp.APBS_MIN_CUTOFF

        ##### REVERSE SIGN OF LOG(ABS(GRID_NEG)) AND DOUBLE BOTH
        logpos *=  2 # this way the range of points varies between
        logneg *= -2 # 2*APBS_MIN_CUTOFF and 2*APBS_MAX_CUTOFF

        ##### RESULT
        self.grid[self.grid > 0] = logpos
        self.grid[self.grid < 0] = logneg

        ##### additional adjustments before saving to json
        self.pack_data()
        self.data["minPotential"] = float(-2 * (vp.APBS_MAX_CUTOFF - vp.APBS_MIN_CUTOFF))
        self.data["maxPotential"] = float( 2 * (vp.APBS_MAX_CUTOFF - vp.APBS_MIN_CUTOFF))
        self.save_data(override_prefix = "apbslog")
        timer.end()


# //////////////////////////////////////////////////////////////////////////////
