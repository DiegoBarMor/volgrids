from ._core.grid_ve import GridVolumetricEnergy

from ._core.app_veins import AppVeins


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

### PLACEHOLDER: default configs go here

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by inherited `AppSubcommand` classes.
### CLI parsed by freyacli.

MODE: str = '' # mode of the application, i.e. "energies"

import pathlib as _pathlib
PATH_STRUCT:       _pathlib.Path = None # "path/input/structure.pdb"
PATH_ENERGIES_CSV: _pathlib.Path = None # "path/input/energies.csv"
PATH_TRAJ:         _pathlib.Path = None # "path/input/traj.xtc"
FOLDER_OUT:        _pathlib.Path = None # "path/output/"

ENERGY_CUTOFF: float # Energies below this cutoff will be ignored (default 1e-3)
