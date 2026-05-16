from ._operations.residues_nucleic import ResiduesNucleic
from ._operations.histogram import Histogram

from ._occupancies.hbaccepts import OgHBAccepts
from ._occupancies.hbdonors import OgHBDonors
from ._occupancies.stacking import OgStacking
from ._occupancies.hydrophobic import OgHydrophobic
# [TODO] hydrophilic? electrostatic?

from ._occupancies.app_occupancy import AppOccupancy
from .app_smutils import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

OG_RADIUS_STACKING: float = 2.0
OG_RADIUS_HBA: float = 2.0
OG_RADIUS_HBD: float = 2.0
OG_RADIUS_HYDROPHOBIC: float = 2.0
# OG_RADIUS_HYDROPHILIC: float = 2.0
# OG_RADIUS_APBS: float = 2.0

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")


######################## UI WITH THE OLD SCHEMA
# ######################## COMMAND LINE ARGUMENTS GLOBALS ########################
# ### These are global variables that are to be set by
# ### an instance of ParamHandler (or its inherited classes)

# import pathlib as _pathlib
# PATH_PDB:    _pathlib.Path = None # "path/input/structure.pdb"
# STR_SMILES:  _pathlib.Path = None # "O=Cc1ccc(O)c(OC)c1"
# PATH_OUTPUT: _pathlib.Path = None # "path/output/table.chem"

# INPUT_KIND: str = '' # either "pdb" or "smiles"
