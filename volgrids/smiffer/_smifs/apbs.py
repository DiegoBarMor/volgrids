import numpy as np

import volgrids as vg
import volgrids.smiffer as sm
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class SmifAPBS(sm.Smif):
    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid) -> vg.Grid:
        if sm.PATH_APBS is not None:
            return self._apbs_to_smif(grid.box, sm.PATH_APBS)

        timer = vg.Timer().start()
        ### "sm.PATH_STRUCT.name" must be used, don't use "self.ms.molname"
        with vg.APBSSubprocess(self.ms.system.atoms, sm.PATH_STRUCT.name) as path_apbs:
            return self._apbs_to_smif(grid.box, path_apbs, timer)


    # --------------------------------------------------------------------------
    @staticmethod
    def _apbs_to_smif(box_dst: vg.Box, path_apbs_in, timer: vg.Timer = None) -> vg.Grid:
        if timer is not None: sm.APBS_ELAPSED_TIME = timer.end(
            text = fy.Color.red("APBS"), end = ' '
        )

        apbs = vg.GridIO.read_auto(path_apbs_in)
        apbs.reshape_as_box(box_dst)

        apbs.dirty = True
        return apbs


    # --------------------------------------------------------------------------
    @staticmethod
    def apply_logabs_transform(grid: vg.Grid) -> None:
        if grid.is_empty():
            print(f"...--- {fy.Color.red('APBS potential grid is empty')}. Skipping logabs transform.", flush = True)
            return

        logpos = np.log10( grid.arr[grid.arr > 0])
        logneg = np.log10(-grid.arr[grid.arr < 0])

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
        grid.arr[grid.arr > 0] = logpos
        grid.arr[grid.arr < 0] = logneg


# //////////////////////////////////////////////////////////////////////////////
