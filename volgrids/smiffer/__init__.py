from .misc.mol_system import MolType, MolSystemSmiffer
from .misc.trimmer import Trimmer

from .parsers.chem_table import ChemTable

from .smifs.smif import Smif
from .smifs.apbs import SmifAPBS
from .smifs.hbonds import SmifHBARing, SmifHBDRing, SmifHBDCone
from .smifs.hydro import SmifHydrophilic, SmifHydrophobic
from .smifs.stacking import SmifStacking

from .ui.param_handler import ParamHandlerSmiffer
from .ui.app import AppSmiffer


############################# CONFIG FILE GLOBALS ##############################
DO_SMIF_STACKING:    bool
DO_SMIF_HBA:         bool
DO_SMIF_HBD:         bool
DO_SMIF_HYDROPHOBIC: bool
DO_SMIF_HYDROPHILIC: bool
DO_SMIF_APBS:        bool

DO_SMIF_LOG_APBS:  bool
DO_SMIF_HYDRODIFF: bool

DO_TRIMMING_SPHERE:    bool
DO_TRIMMING_OCCUPANCY: bool
DO_TRIMMING_RNDS:      bool
DO_TRIMMING_FARAWAY:   bool
SAVE_TRIMMING_MASK:    bool

TRIMMING_DIST_SMALL: float
TRIMMING_DIST_MID:   float
TRIMMING_DIST_LARGE: float

MAX_RNDS_DIST:   float
COG_CUBE_RADIUS: int

TRIM_FARAWAY_DIST: float

APBS_MIN_CUTOFF: int
APBS_MAX_CUTOFF: int

ENERGY_SCALE: float

MU_HYDROPHOBIC:    float
SIGMA_HYDROPHOBIC: float

MU_HYDROPHILIC:    float
SIGMA_HYDROPHILIC: float

MU_ANGLE_HBA:    float
MU_DIST_HBA:     float
SIGMA_ANGLE_HBA: float
SIGMA_DIST_HBA:  float

MU_ANGLE_HBD:    float
MU_DIST_HBD:     float
SIGMA_ANGLE_HBD: float
SIGMA_DIST_HBD:  float

MU_ANGLE_HBD_FIXED:    float
MU_DIST_HBD_FIXED:     float
SIGMA_ANGLE_HBD_FIXED: float
SIGMA_DIST_HBD_FIXED:  float

MU_ANGLE_STACKING: float
MU_DIST_STACKING:  float
COV_STACKING_00:   float
COV_STACKING_01:   float
COV_STACKING_10:   float
COV_STACKING_11:   float

GAUSSIAN_KERNEL_SIGMAS: int

__config_keys__ = set(__annotations__.keys())


############################### NUMERIC GLOBALS ################################
import numpy as _np
MU_HBA:      _np.ndarray[float, float]
COV_HBA:     _np.ndarray
COV_INV_HBA: _np.ndarray

MU_HBD:      _np.ndarray[float, float]
COV_HBD:     _np.ndarray
COV_INV_HBD: _np.ndarray

MU_HBD_FIXED:      _np.ndarray[float, float]
COV_HBD_FIXED:     _np.ndarray
COV_INV_HBD_FIXED: _np.ndarray

MU_STACKING:      _np.ndarray[float, float]
COV_STACKING:     _np.ndarray
COV_INV_STACKING: _np.ndarray

SIGMA_DIST_STACKING: float


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

import pathlib as _pathlib
PATH_STRUCTURE:  _pathlib.Path = None # "path/input/struct.pdb"
PATH_TRAJECTORY: _pathlib.Path = None # "path/input/traj.xtc"
PATH_APBS:       _pathlib.Path = None # "path/input/apbs.pqr.dx"
PATH_TABLE:      _pathlib.Path = None # "path/input/table.chem"
FOLDER_OUT:      _pathlib.Path = None # "folder/output/"

PS_INFO: tuple[float, float, float, float] = None # pocket sphere info: [radius, x, y, z]
CURRENT_MOLTYPE: MolType = MolType.NONE           # type of the current molecule

USE_STRUCTURE_HYDROGENS = False # whether to use hydrogens from the structure to calculate hbond smifs
