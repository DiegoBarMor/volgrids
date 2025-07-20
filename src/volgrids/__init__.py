from ._framework._core.grid import Grid
from ._framework._core.mol_system import MolSystem

from ._framework._kernels.kernel import Kernel
from ._framework._kernels.boolean import \
    KernelSphere, KernelCylinder, KernelDisk, KernelDiskConecut
from ._framework._kernels.gaussian import \
    KernelGaussianUnivariateDist, KernelGaussianBivariateAngleDist

from ._framework._misc.math import Math
from ._framework._misc.params_gaussian import \
    ParamsGaussianUnivariate, ParamsGaussianBivariate
from ._framework._misc.timer import Timer
from ._framework._misc.utils import resolve_path

from ._framework._parsers.parser_ini import ParserIni
from ._framework._parsers.parser_config import ParserConfig
from ._framework._parsers.grid_io import GridFormat, GridIO

from ._framework._ui.param_handler import ParamHandler
from ._framework._ui.app import App


############################# CONFIG FILE GLOBALS ##############################
OUTPUT_FORMAT: GridFormat

GZIP_COMPRESSION: int
FLOAT_DTYPE: type
WARNING_GRID_SIZE: float

GRID_DX: float
GRID_DY: float
GRID_DZ: float

GRID_XRES: int
GRID_YRES: int
GRID_ZRES: int

EXTRA_BOX_SIZE: int
USE_FIXED_DELTAS: bool

__config_keys__ = set(__annotations__.keys())


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

import pathlib as _pathlib
PATH_CUSTOM_CONFIG: _pathlib.Path = None # "path/input/globals.ini"
