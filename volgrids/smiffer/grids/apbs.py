import numpy as np
import volgrids.vgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class GridAPBS(vg.Grid):
    def get_type(self):
        return "apbs"

    def run(self, trimmer: "sm.GridTrimmer"):
        if self.ms.meta.path_apbs is None:
            print("...--- APBS output not provided. Electrostatic potential skipped.", flush = True)
            return
        super().run(trimmer)

    def populate_grid(self):
        xmin1, ymin1, zmin1 = self.ms.minCoords
        xmax1, ymax1, zmax1 = self.ms.maxCoords
        xres1, yres1, zres1 = self.ms.resolution
        _,apbs = vg.read_auto(self.ms.meta.path_apbs)

        self.grid = vg.interpolate_3d(
            x0 = np.linspace(apbs.xmin, apbs.xmax, apbs.xres),
            y0 = np.linspace(apbs.ymin, apbs.ymax, apbs.yres),
            z0 = np.linspace(apbs.zmin, apbs.zmax, apbs.zres),
            data_0 = apbs.grid,
            new_coords = np.mgrid[
                xmin1 : xmax1 : complex(0, xres1),
                ymin1 : ymax1 : complex(0, yres1),
                zmin1 : zmax1 : complex(0, zres1),
            ].T
        ).astype(vg.FLOAT_DTYPE)


    def apply_logabs_transform(self):
        if self.is_empty():
            print(f"...--- APBS potential grid is empty. Skipping logabs transform.", flush = True)
            return

        logpos = np.log10( self.grid[self.grid > 0])
        logneg = np.log10(-self.grid[self.grid < 0])

        ##### APPLY CUTOFFS
        logpos[logpos < sm.APBS_MIN_CUTOFF] = sm.APBS_MIN_CUTOFF
        logneg[logneg < sm.APBS_MIN_CUTOFF] = sm.APBS_MIN_CUTOFF
        logpos[logpos > sm.APBS_MAX_CUTOFF] = sm.APBS_MAX_CUTOFF
        logneg[logneg > sm.APBS_MAX_CUTOFF] = sm.APBS_MAX_CUTOFF

        ##### SHIFT VALUES TO 0
        logpos -= sm.APBS_MIN_CUTOFF
        logneg -= sm.APBS_MIN_CUTOFF

        ##### REVERSE SIGN OF LOG(ABS(GRID_NEG)) AND DOUBLE BOTH
        logpos *=  2 # this way the range of points varies between
        logneg *= -2 # 2*APBS_MIN_CUTOFF and 2*APBS_MAX_CUTOFF

        ##### RESULT
        self.grid[self.grid > 0] = logpos
        self.grid[self.grid < 0] = logneg
        self.save_data(override_prefix = "apbslog")


# //////////////////////////////////////////////////////////////////////////////
