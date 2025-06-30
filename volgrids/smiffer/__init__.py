from .grids.apbs import GridAPBS
# from .grids.cavities import GridCavities
from .grids.hbonds import GridHBAccepts, GridHBDonors
from .grids.hydro import GridHydrophilic, GridHydrophobic
from .grids.stacking import GridStacking

from .parsers.args_parser import SmifferArgsParser
from .parsers.chem_table import ChemTable

from .trimmers import GridTrimmer
from .systems import MolType, SmifferMolecularSystem

from .smiffer import SmifferCalculator


############################# CONFIG FILE GLOBALS ##############################
DO_SMIF_STACKING: bool
DO_SMIF_HBA: bool
DO_SMIF_HBD: bool
DO_SMIF_HYDROPHOBIC: bool
DO_SMIF_HYDROPHILIC: bool
DO_SMIF_APBS: bool

DO_SMIF_LOG_APBS: bool
DO_SMIF_HYDRODIFF: bool

DO_TRIMMING_SPHERE: bool
DO_TRIMMING_OCCUPANCY: bool
DO_TRIMMING_RNDS: bool
SAVE_CACHED_MASK: bool

TRIMMING_DIST_LARGE: float
TRIMMING_DIST_SMALL: float

MAX_RNDS_DIST: float
COG_CUBE_RADIUS: int

APBS_MIN_CUTOFF: int
APBS_MAX_CUTOFF: int

MU_HYDROPHOBIC: float
SIGMA_HYDROPHOBIC: float

MU_HYDROPHILIC: float
SIGMA_HYDROPHILIC: float

MU_ANGLE_HBA: float
MU_ANGLE_HBD: float
MU_DIST_HBOND: float

SIGMA_DIST_HBOND: float
SIGMA_ANGLE_HBA: float
SIGMA_ANGLE_HBD: float

MU_ANGLE_STACKING: float
MU_DIST_STACKING: float

COV_STACKING_00: float
COV_STACKING_01: float
COV_STACKING_10: float
COV_STACKING_11: float

GAUSSIAN_KERNEL_SIGMAS: int


import numpy as _np
import volgrids as _vg

_vg.ConfigParser(_vg.resolve_path("volgrids/config.ini")).apply_config(
    key = "SMIFFER",
    globs = globals(),
    valid_configs = set(__annotations__.keys())
)


############################### NUMERIC GLOBALS ################################
MU_HBA = _np.array([MU_ANGLE_HBA, MU_DIST_HBOND])
COV_HBA = _np.array(
    [[SIGMA_ANGLE_HBA**2, 0],
    [0, SIGMA_DIST_HBOND**2]]
)
COV_INV_HBA = _np.linalg.inv(COV_HBA)


MU_HBD = _np.array([MU_ANGLE_HBD, MU_DIST_HBOND])
COV_HBD =  _np.array(
    [[SIGMA_ANGLE_HBD**2, 0],
    [0, SIGMA_DIST_HBOND**2]]
)
COV_INV_HBD = _np.linalg.inv(COV_HBD)


MU_STACKING = _np.array([MU_ANGLE_STACKING, MU_DIST_STACKING])
COV_STACKING = _np.array(
    [[COV_STACKING_00, COV_STACKING_01],
     [COV_STACKING_10, COV_STACKING_11]]
)
COV_INV_STACKING = _np.linalg.inv(COV_STACKING)


### square root of the DIST contribution to COV_STACKING,
SIGMA_DIST_STACKING = _np.sqrt(COV_STACKING_11)


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ArgsParser (or its inherited classes)

from pathlib import Path
PATH_STRUCTURE:  Path = None # "path/input/struct.pdb"
PATH_TRAJECTORY: Path = None # "path/input/traj.xtc"
PATH_APBS:       Path = None # "path/input/apbs.pqr.dx"
PATH_TABLE:      Path = None # "path/input/table.chem"
PATH_CONFIG:     Path = None # "path/input/globals.ini"
FOLDER_OUT:      Path = None # "folder/output/"

PS_INFO: tuple[float, float, float, float] = None # pocket sphere info: [radius, x, y, z]
CURRENT_MOLTYPE: MolType = MolType.NONE           # type of the current molecule
