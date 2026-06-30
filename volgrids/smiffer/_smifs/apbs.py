import numpy as np

import volgrids as vg
import volgrids.smiffer as smf
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class SmifAPBS(smf.Smif):
    # --------------------------------------------------------------------------
    def populate_grid(self, grid: vg.Grid) -> None:
        if smf.PATH_APBS is not None:
            return self._apbs_to_smif(grid, smf.PATH_APBS)

        timer = vg.Timer().start()
        ### "smf.PATH_STRUCT.name" must be used, don't use "self.mm.molname"
        with vg.APBSSubprocess(self.mm.get_all_atoms(), smf.PATH_STRUCT.name) as path_apbs:
            return self._apbs_to_smif(grid, path_apbs, timer)


    # --------------------------------------------------------------------------
    def gen_pqr(self):
        with vg.APBSSubprocess(self.mm.get_all_atoms(), smf.PATH_STRUCT.name, only_pdb2pqr = True) as _:
            return


    # --------------------------------------------------------------------------
    @staticmethod
    def _apbs_to_smif(grid: vg.Grid, path_apbs_in, timer: vg.Timer = None) -> None:
        if timer is not None: smf.APBS_ELAPSED_TIME = timer.end(
            text = fy.Color.red("APBS"), end = ' '
        )

        apbs = vg.Grid.load(path_apbs_in)
        apbs.reshape_as_box(grid.box)
        grid.arr = apbs.arr
        grid.dirty = True
        del apbs


    # --------------------------------------------------------------------------
    @staticmethod
    def apply_logabs_transform(grid: vg.Grid) -> None:
        if grid.is_empty():
            print(f"...--- {fy.Color.red('APBS potential grid is empty')}. Skipping logabs transform.", flush = True)
            return

        logpos = np.log10( grid.arr[grid.arr > 0])
        logneg = np.log10(-grid.arr[grid.arr < 0])

        ##### APPLY CUTOFFS
        logpos[logpos < vg.CFG.misc_logapbs_min] = vg.CFG.misc_logapbs_min
        logneg[logneg < vg.CFG.misc_logapbs_min] = vg.CFG.misc_logapbs_min
        logpos[logpos > vg.CFG.misc_logapbs_max] = vg.CFG.misc_logapbs_max
        logneg[logneg > vg.CFG.misc_logapbs_max] = vg.CFG.misc_logapbs_max

        ##### SHIFT VALUES TO 0
        logpos -= vg.CFG.misc_logapbs_min
        logneg -= vg.CFG.misc_logapbs_min

        ##### REVERSE SIGN OF LOG(ABS(GRID_NEG)) AND DOUBLE BOTH
        logpos *=  2 # this way the range of points varies between
        logneg *= -2 # 2*MISC_LOGAPBS_MIN and 2*MISC_LOGAPBS_MAX

        ##### RESULT
        grid.arr[grid.arr > 0] = logpos
        grid.arr[grid.arr < 0] = logneg


# //////////////////////////////////////////////////////////////////////////////
