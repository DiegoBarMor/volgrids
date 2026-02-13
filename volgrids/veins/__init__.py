from ._grids.grid_ve import GridVolumetricEnergy
from ._grids.grid_vf import GridVolumetricForce

from ._ui.param_handler import ParamHandlerVeins
from ._ui.app import AppVeins


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

ENERGY_CUTOFF:    float = 1e-3 # Energies below this cutoff will be ignored
FORCE_CUTOFF:     float = 1e-3 # Force magnitudes below this cutoff will be ignored
FORCE_VECTOR_LEN: float = 5.0  # Length of the strongest force vectors when visualized (in Angstroms)

__config_keys__ = set(globals().keys()) - _keys_other; del _keys_other


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

MODE: str = '' # mode of the application, i.e. "energies", "forces"
TRAJ_FRAME: int | None = None # current trajectory frame being processed (only relevant in trajectory mode)

import pathlib as _pathlib
PATH_CSV_IN: _pathlib.Path | None = None # "path/input/energies.csv"
FOLDER_OUT:  _pathlib.Path | None = None # "path/output/"
