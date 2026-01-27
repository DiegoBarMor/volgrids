from ._core.mol_system import MolType, MolSystemSmiffer
from ._core.trimmer import Trimmer

from ._parsers.parser_chem_table import ParserChemTable

from ._smifs.smif import Smif
from ._smifs.apbs import SmifAPBS
from ._smifs._hbonds.hbaccepts import SmifHBAccepts
from ._smifs._hbonds.hbdonors import SmifHBDonors
from ._smifs._hydro.hydrophilic import SmifHydrophilic
from ._smifs._hydro.hydrophobic import SmifHydrophobic
from ._smifs.stacking import SmifStacking

from ._ui.param_handler import ParamHandlerSmiffer
from ._ui.app import AppSmiffer


############################# CONFIG FILE GLOBALS ##############################
import numpy as _np
_keys_other = set(globals().keys())

DO_SMIF_STACKING:    bool = True
DO_SMIF_HBA:         bool = True
DO_SMIF_HBD:         bool = True
DO_SMIF_HYDROPHOBIC: bool = True
DO_SMIF_HYDROPHILIC: bool = True
DO_SMIF_APBS:        bool = True

DO_SMIF_LOG_APBS:  bool = False
DO_SMIF_HYDRODIFF: bool = False

DO_TRIMMING_SPHERE:    bool = True
DO_TRIMMING_OCCUPANCY: bool = True
DO_TRIMMING_RNDS:      bool = True
DO_TRIMMING_FARAWAY:   bool = True
SAVE_TRIMMING_MASK:    bool = False

USE_STRUCTURE_HYDROGENS: bool = True

TRIMMING_DIST_SMALL: float = 2.5
TRIMMING_DIST_MID:   float = 3.0
TRIMMING_DIST_LARGE: float = 3.5

MAX_RNDS_DIST:   float = _np.inf
COG_CUBE_RADIUS: int = 4

TRIM_FARAWAY_DIST: float = 7.0

ENERGY_SCALE: float = 3.5

MU_HYDROPHOBIC:    float = 4.4
SIGMA_HYDROPHOBIC: float = 0.3

MU_HYDROPHILIC:    float = 3.0
SIGMA_HYDROPHILIC: float = 0.15

MU_ANGLE_HBA:    float = 129.9
MU_DIST_HBA:     float = 2.93
SIGMA_ANGLE_HBA: float = 20.0
SIGMA_DIST_HBA:  float = 0.21

MU_ANGLE_HBD_FREE:    float = 109.0
MU_DIST_HBD_FREE:     float = 2.93
SIGMA_ANGLE_HBD_FREE: float = 20.0
SIGMA_DIST_HBD_FREE:  float = 0.21

MU_ANGLE_HBD_FIXED:    float = 180.0
MU_DIST_HBD_FIXED:     float = 2.93
SIGMA_ANGLE_HBD_FIXED: float = 30.0
SIGMA_DIST_HBD_FIXED:  float = 0.21

MU_ANGLE_STACKING: float = 29.97767535
MU_DIST_STACKING:  float = 4.1876158
COV_STACKING_00:   float = 169.9862228
COV_STACKING_01:   float = 6.62318852
COV_STACKING_10:   float = 6.62318852
COV_STACKING_11:   float = 0.37123882


GAUSSIAN_KERNEL_SIGMAS: int = 4
APBS_MIN_CUTOFF: int = -2
APBS_MAX_CUTOFF: int = 3

__config_keys__ = set(globals().keys()) - _keys_other


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

import pathlib as _pathlib
PATH_STRUCT: _pathlib.Path = None # "path/input/struct.pdb"
PATH_TRAJ:   _pathlib.Path = None # "path/input/traj.xtc"
PATH_APBS:   _pathlib.Path = None # "path/input/apbs.pqr.dx"
PATH_TABLE:  _pathlib.Path = None # "path/input/table.chem"
FOLDER_OUT:  _pathlib.Path = None # "folder/output/"

SPHERE_INFO: tuple[float, float, float, float] = None # pocket sphere info: [x, y, z, radius]
CURRENT_MOLTYPE: MolType = MolType.NONE               # type of the current molecule


############################### RUNTIME GLOBALS ################################
import volgrids as _vg
PARAMS_HBA:       _vg.ParamsGaussianBivariate
PARAMS_HBD_FREE:  _vg.ParamsGaussianBivariate
PARAMS_HBD_FIXED: _vg.ParamsGaussianBivariate
PARAMS_HPHOB:     _vg.ParamsGaussianUnivariate
PARAMS_HPHIL:     _vg.ParamsGaussianUnivariate
PARAMS_STACK:     _vg.ParamsGaussianBivariate
SIGMA_DIST_STACKING: float

APBS_ELAPSED_TIME: float = 0.0
