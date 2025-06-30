from .args_parser import VeinsArgsParser
from .grid_ve import GridVolumetricEnergy
from .veins import VeinsApp

######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ArgsParser (or its inherited classes)

from pathlib import Path
PATH_STRUCTURE:    Path = None  # "path/input/structure.pdb"
PATH_ENERGIES_CSV: Path = None  # "path/input/energies.csv"
FOLDER_OUT:        Path = None  # "path/output/"

ENERGY_CUTOFF:     float = 1e-3  # Energies below this cutoff will be ignored
