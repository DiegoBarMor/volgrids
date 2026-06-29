from ._core.mol_system import MolSystem
from ._core.traj_multiprocess import TrajMultiprocess
from ._core.cavity_finder import CavityFinder
from ._core.trimmer import Trimmer

from ._parsers.resname_standard import ResnameStandard
from ._parsers.parser_chem_table import ParserChemTable

from ._smifs import _core as _smifs_core
from ._smifs._core import Smif
from ._smifs.apbs import SmifAPBS
from ._smifs.hbaccepts import SmifHBAccepts
from ._smifs.hbdonors import SmifHBDonors
from ._smifs.hydrophilic import SmifHydrophilic
from ._smifs.hydrophobic import SmifHydrophobic
from ._smifs.stacking import SmifStacking

from .app_smiffer import AppSmiffer


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

SMIF_STK:   bool = True
SMIF_HBA:   bool = True
SMIF_HBD:   bool = True
SMIF_APBS:  bool = True
SMIF_HPHOB: bool = True
SMIF_HPHIL: bool = True
SMIF_USE_HYDROGENS: bool = True
SMIF_HB_ONLY_NBASE: bool = False

TRIM_SPHERE:    bool = True
TRIM_OCCUPANCY: bool = True
TRIM_RNDS:      bool = True
TRIM_FARAWAY:   bool = True
TRIM_CAVITIES:  bool = False
TRIM_SAVE: bool = False

CAV_THRESHOLD: int = 3
CAV_NPASSES: int = 2
CAV_WEIGHT: float = 0.0
CAV_SAVE: bool = False

TRIM_OCC_DIST_SHORT: float = 2.5
TRIM_OCC_DIST_MID:   float = 3.0
TRIM_OCC_DIST_LONG:  float = 3.5

TRIM_RNDS_MAX_DIST:    float = float("inf")
TRIM_RNDS_CUBE_RADIUS: int = 4

TRIM_FARAWAY_DIST: float = 7.0

PARAM_STK_SCALE: float = 3.5
PARAM_HB_SCALE:  float = 3.5

PARAM_HPHOB_DIST_MU:    float = 4.4
PARAM_HBHOB_DIST_SIGMA: float = 0.3

PARAM_HPHIL_DIST_MU:    float = 3.0
PARAM_HPHIL_DIST_SIGMA: float = 0.15

PARAM_HBA_ANGLE_MU:    float = 129.9
PARAM_HBA_DIST_MU:     float = 2.93
PARAM_HBA_ANGLE_SIGMA: float = 20.0
PARAM_HBA_DIST_SIGMA:  float = 0.21

PARAM_HBD_FREE_ANGLE_MU:    float = 109.0
PARAM_HBD_FREE_DIST_MU:     float = 2.93
PARAM_HBD_FREE_ANGLE_SIGMA: float = 20.0
PARAM_HBD_FREE_DIST_SIGMA:  float = 0.21

PARAM_HBD_FIXED_ANGLE_MU:    float = 180.0
PARAM_HBD_FIXED_DIST_MU:     float = 2.93
PARAM_HBD_FIXED_ANGLE_SIGMA: float = 30.0
PARAM_HBD_FIXED_DIST_SIGMA:  float = 0.21

PARAM_STK_ANGLE_MU: float = 29.97767535
PARAM_STK_DIST_MU:  float = 4.1876158
PARAM_STK_COV00:    float = 169.9862228
PARAM_STK_COV01:    float = 6.62318852
PARAM_STK_COV10:    float = 6.62318852
PARAM_STK_COV11:    float = 0.37123882

MISC_KERNEL_GAUSSIAN_SIGMAS: int = 4
MISC_LOGAPBS_MIN: int = -2
MISC_LOGAPBS_MAX: int = 3

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by inherited `AppSubcommand` classes.
### CLI parsed by freyacli.

import pathlib as _pathlib
import volgrids as _vg
PATH_STRUCT:      _pathlib.Path = None # "path/input/struct.pdb"
PATH_APBS:        _pathlib.Path = None # "path/input/apbs.pqr.dx"
PATH_CHEM_LIGAND: _pathlib.Path = None # "path/input/table.chem"

SPHERES: list[_vg.SphereInfo] = [] # list of pocket sphere infos: [[x, y, z, radius], ...]
BOXES_ENFORCED: list[_vg.Box] = [] # list of boxes enforced by the user: [[x_min, x_max, y_min, y_max, z_min, z_max], ...]

CUSTOM_RESIDUES: str = "" # "A.3 A.4 A.5 B.10 ..."


############################### RUNTIME GLOBALS ################################
PARAMS_HBA:       _vg.ParamsGaussianBivariate
PARAMS_HBD_FREE:  _vg.ParamsGaussianBivariate
PARAMS_HBD_FIXED: _vg.ParamsGaussianBivariate
PARAMS_HPHOB:     _vg.ParamsGaussianUnivariate
PARAMS_HPHIL:     _vg.ParamsGaussianUnivariate
PARAMS_STACK:     _vg.ParamsGaussianBivariate
SIGMA_DIST_STACKING: float

APBS_ELAPSED_TIME: float = 0.0
