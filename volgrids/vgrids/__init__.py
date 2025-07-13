from .kernels.kernel import Kernel
from .kernels.boolean import \
    KernelSphere, KernelCylinder, KernelDisk, KernelDiskConecut
from .kernels.gaussian import \
    KernelGaussianUnivariate, KernelGaussianMultivariate

from .misc.grid import Grid
from .misc.math import Math
from .misc.mol_system import MolSystem
from .misc.utils import Timer, resolve_path

from .parsers.parser_ini import ParserIni
from .parsers.parser_config import ParserConfig
from .parsers.grid_io import GridFormat, GridIO

from .ui.param_handler import ParamHandler
from .ui.app import App


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


import numpy as np  # [TODO] move these lines
ParserConfig(resolve_path("volgrids/config.ini")).apply_config(
    key = "VGRIDS", scope = globals(),
    valid_configs = set(__annotations__.keys())
)
