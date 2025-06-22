import numpy as _np

### OUTPUTS
DO_OUTPUT_DX = False
DO_OUTPUT_MRC = True
DO_OUTPUT_CMAP = False

######################## SPACE EFFICIENCY
GZIP_COMPRESSION = 9 # gzip compression level for CMAP files (0-9); h5py default is 4
FLOAT_DTYPE = _np.float32 # numerical precision of the grid data
WARNING_GRID_SIZE = 5.0e7 # if the grid would exceed this amount of points, trigger a warning with possibility to abort

######################## GRIDS
### deltas used for calculations when USE_FIXED_DELTAS = True (resolutions change)
GRID_DX = 0.25
GRID_DY = 0.25
GRID_DZ = 0.25

### resolution used for calculations when USE_FIXED_DELTAS = False (deltas change)
GRID_XRES = 200
GRID_YRES = 200
GRID_ZRES = 200

EXTRA_BOX_SIZE = 5 # only applies to whole mode
USE_FIXED_DELTAS = True # whether to use fixed dx,dy,dz and le xres,yres,zres change (or the opposite)


#################################################################

from ._core.args import ArgsParser
from ._core.grid import Grid
from ._core.kernel import Kernel
from ._core.systems import MolecularSystem

from ._core.kernels.boolean import \
    KernelSphere, KernelCylinder, KernelDisk, KernelDiskConecut
from ._core.kernels.gaussian import \
    KernelGaussianUnivariate, KernelGaussianMultivariate

from ._core.utils.math import \
    normalize, dot_product, get_norm, get_angle, \
    get_projection, get_projection_height, \
    univariate_gaussian, multivariate_gaussian, \
    interpolate_3d, format_vector_str, get_coords_array
from ._core.utils.ini_parser import IniParser
from ._core.utils.io import \
    read_dx, read_mrc, read_cmap, \
    write_dx, write_mrc, write_cmap, \
    read_auto, get_cmap_keys, resolve_path
from ._core.utils.timer import Timer
