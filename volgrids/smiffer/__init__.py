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
