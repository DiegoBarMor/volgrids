from ._core.mol_type import MolType
from ._core.mol_system import MolSystemSmiffer
from ._core.cavity_finder import CavityFinder
from ._core.trimmer import Trimmer

from ._parsers.parser_chem_table import ParserChemTable

from ._smifs.smif import Smif
from ._smifs.apbs import SmifAPBS
from ._smifs._hbonds.hbaccepts import SmifHBAccepts
from ._smifs._hbonds.hbdonors import SmifHBDonors
from ._smifs._hydro.hydrophilic import SmifHydrophilic
from ._smifs._hydro.hydrophobic import SmifHydrophobic
from ._smifs.stacking import SmifStacking

# Probe-Probe (PP) fields - spherical accessibility regions
from ._smifs._hbonds.hbaccepts_pp import SmifHBAcceptsPP
from ._smifs._hbonds.hbdonors_pp import SmifHBDonorsPP
from ._smifs.stacking_pp import SmifStackingPP
from ._smifs._hydro.hydro_pp import SmifHydroPP

from ._misc.sphere_info import SphereInfo

from ._ui.param_handler import ParamHandlerSmiffer
from ._ui.app import AppSmiffer


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

GRID_FORMAT_OUTPUT: str = "CMAP_PACKED"

DO_SMIF_STACKING:    bool = True
DO_SMIF_HBA:         bool = True
DO_SMIF_HBD:         bool = True
DO_SMIF_HYDROPHOBIC: bool = True
DO_SMIF_HYDROPHILIC: bool = True
DO_SMIF_APBS:        bool = True

# Probe-Probe (PP) fields - spherical accessibility regions
DO_SMIF_PP:          bool = False  # Generate all PP fields
DO_SMIF_HBA_PP:      bool = False  # HB acceptor probe-probe (hbaPP)
DO_SMIF_HBD_PP:      bool = False  # HB donor probe-probe (hbdPP)
DO_SMIF_STACKING_PP: bool = False  # Stacking probe-probe (stkPP)
DO_SMIF_HYDRO_PP:    bool = False  # Hydrophobic probe-probe (hpPP)

DO_SMIF_LOG_APBS:  bool = False
DO_SMIF_HYDRODIFF: bool = False

DO_TRIMMING_SPHERE:    bool = True
DO_TRIMMING_OCCUPANCY: bool = True
DO_TRIMMING_RNDS:      bool = True
DO_TRIMMING_FARAWAY:   bool = True
DO_TRIMMING_CAVITIES:  bool = False

TRIMMING_CAVITIES_THRESHOLD: int = 3
CAVITIES_NPASSES: int = 2
CAVITIES_WEIGHT: float = 0.0

SAVE_TRIMMING_MASK: bool = False
SAVE_CAVITIES: bool = False

USE_STRUCTURE_HYDROGENS: bool = True
HBONDS_ONLY_NUCLEOBASE: bool = False

TRIMMING_DIST_TINY:  float = 1.0
TRIMMING_DIST_SMALL: float = 2.5
TRIMMING_DIST_MID:   float = 3.0
TRIMMING_DIST_LARGE: float = 3.5

MAX_RNDS_DIST:   float = float("inf")
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
__config_keys__.remove("_keys_other")


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

import pathlib as _pathlib
PATH_STRUCT: _pathlib.Path = None # "path/input/struct.pdb"
PATH_TRAJ:   _pathlib.Path = None # "path/input/traj.xtc"
PATH_APBS:   _pathlib.Path = None # "path/input/apbs.pqr.dx"
PATH_TABLE:  _pathlib.Path = None # "path/input/table.chem"
FOLDER_OUT:  _pathlib.Path = None # "folder/output/"

SPHERE: SphereInfo  = None              # pocket sphere info: [x, y, z, radius]
CURRENT_MOLTYPE: MolType = MolType.NONE # type of the current molecule

CUSTOM_RESIDS: str = "" # "1 2 4 5 ..."


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
