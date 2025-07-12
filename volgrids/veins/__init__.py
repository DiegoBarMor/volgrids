from .misc.grid_ve import GridVolumetricEnergy

from .ui.param_handler import ParamHandlerVeins
from .ui.app import AppVeins


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

MODE: str = '' # mode of the application, i.e. "energies"

from pathlib import Path
PATH_STRUCTURE:    Path = None # "path/input/structure.pdb"
PATH_ENERGIES_CSV: Path = None # "path/input/energies.csv"
PATH_TRAJECTORY:   Path = None # "path/input/traj.xtc"
FOLDER_OUT:        Path = None # "path/output/"

ENERGY_CUTOFF:     float = 1e-3  # Energies below this cutoff will be ignored
