from abc import ABC, abstractmethod

import volgrids.vgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Smif(vg.Grid, ABC):
    # --------------------------------------------------------------------------
    @abstractmethod
    def populate_grid(self):
        return


# //////////////////////////////////////////////////////////////////////////////
