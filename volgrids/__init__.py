from ._core.grid import Grid
from ._core.mol_system import MolSystem

from ._kernels.kernel import Kernel
from ._kernels.boolean import \
    KernelSphere, KernelCylinder, KernelDisk, KernelDiskConecut
from ._kernels.gaussian import \
    KernelGaussianUnivariateDist, KernelGaussianBivariateAngleDist

from ._misc.math import Math
from ._misc.params_gaussian import ParamsGaussian, \
    ParamsGaussianUnivariate, ParamsGaussianBivariate
from ._misc.timer import Timer
from ._misc.utils import resolve_path_package, resolve_path_resource

from ._parsers.parser_ini import ParserIni
from ._parsers.parser_config import ParserConfig
from ._parsers.grid_io import GridFormat, GridIO

from ._ui.param_handler import ParamHandler
from ._ui.app import App


############################# CONFIG FILE GLOBALS ##############################
import numpy as _np
OUTPUT_FORMAT: GridFormat = GridFormat.CMAP_PACKED

GZIP_COMPRESSION: int = 9
FLOAT_DTYPE: type = _np.float32
WARNING_GRID_SIZE: float = 5.0e7

GRID_DX: float = 0.25
GRID_DY: float = 0.25
GRID_DZ: float = 0.25

GRID_XRES: int = 200
GRID_YRES: int = 200
GRID_ZRES: int = 200

EXTRA_BOX_SIZE: int = 5
USE_FIXED_DELTAS: bool = True

__config_keys__ = set(__annotations__.keys())


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by reading config files
### DEFAULT config.ini allows to first read "config_volgrids.ini" from the volgrid's repo root,
### to be used by the volgrid's main scripts. Its default remains None for any other use case.
### CUSTOM config.ini allows the user to specify a custom config file path from the command line.

import pathlib as _pathlib
PATH_DEFAULT_CONFIG: _pathlib.Path = None # "./config_volgrids.ini"
PATH_CUSTOM_CONFIG:  _pathlib.Path = None # "path/input/config.ini"
