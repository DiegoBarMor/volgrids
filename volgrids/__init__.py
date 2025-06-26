from ._core.grid import Grid
from ._core.kernel import Kernel
from ._core.systems import MolecularSystem
from ._core.math import Math
from ._core.utils import Timer, \
    format_vector_str, resolve_path

from ._core.kernels.boolean import \
    KernelSphere, KernelCylinder, KernelDisk, KernelDiskConecut
from ._core.kernels.gaussian import \
    KernelGaussianUnivariate, KernelGaussianMultivariate

from ._core.parsers.args_parser import ArgsParser
from ._core.parsers.ini_parser import IniParser
from ._core.parsers.config_parser import ConfigParser
from ._core.parsers.grid_io import GridFormat, GridIO


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


ConfigParser(resolve_path("vgrids.config")).apply_config(
    key = "VOLGRIDS",
    globs = globals(),
    valid_configs = set(__annotations__.keys())
)

OUTPUT_FORMAT = GridFormat.from_str(OUTPUT_FORMAT) # [TODO] improve IniParser.parse_str to avoid this
