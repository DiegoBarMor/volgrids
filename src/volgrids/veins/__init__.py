from ._core.grid_ve import GridVolumetricEnergy

from ._ui.param_handler import ParamHandlerVeins
from ._ui.app import AppVeins

__config_keys__ = set(__annotations__.keys())


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

MODE: str = '' # mode of the application, i.e. "energies"

import pathlib as _pathlib
PATH_STRUCTURE:    _pathlib.Path = None # "path/input/structure.pdb"
PATH_ENERGIES_CSV: _pathlib.Path = None # "path/input/energies.csv"
PATH_TRAJECTORY:   _pathlib.Path = None # "path/input/traj.xtc"
FOLDER_OUT:        _pathlib.Path = None # "path/output/"

ENERGY_CUTOFF: float # Energies below this cutoff will be ignored (default 1e-3)
