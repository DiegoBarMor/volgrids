import numpy as np
import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class GridAPBS(vg.Grid):
    # --------------------------------------------------------------------------
    def get_type(self):
        return "apbs"


    # --------------------------------------------------------------------------
    def populate_grid(self):
        _,apbs = vg.GridIO.read_auto(self.ms.meta.path_apbs)
        apbs.reshape(
            new_min = (self.xmin, self.ymin, self.zmin),
            new_max = (self.xmax, self.ymax, self.zmax),
            new_res = (self.xres, self.yres, self.zres)
        )
        self.grid = apbs.grid


    # --------------------------------------------------------------------------
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


# //////////////////////////////////////////////////////////////////////////////
