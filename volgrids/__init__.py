################################################################################
############################### GLOBAL VARIABLES ###############################
################################################################################
import numpy as _np

######################## TOGGLES
### GRIDS
DO_SMIF_STACKING = True
DO_SMIF_HBA = True
DO_SMIF_HBD = True
DO_SMIF_HYDROPHOBIC = True
DO_SMIF_HYDROPHILIC = True
DO_SMIF_APBS = True

DO_SMIF_LOG_APBS = False
DO_SMIF_HYDRODIFF = False

### TRIMMING
DO_TRIMMING_SPHERE = True
DO_TRIMMING_OCCUPANCY = True
DO_TRIMMING_RNDS = False

### OUTPUTS
DO_OUTPUT_DX = False
DO_OUTPUT_MRC = True
DO_OUTPUT_CMAP = False
SAVE_CACHED_MASK = False # saves the logical inverse of the trimming mask

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

######################## TRIMMING
### OCCUPANCY TRIMMING
TRIMMING_DIST_LARGE = 3.0
TRIMMING_DIST_SMALL = 2.5

### RANDOM SEARCH TRIMMING
MAX_RNDS_DIST = _np.inf
COG_CUBE_RADIUS = 3

######################## LOG_ABS(APBS)
APBS_MIN_CUTOFF = -2
APBS_MAX_CUTOFF =  3

######################## PARAMETERS FOR THE STATISTICAL ENERGY FUNCTIONS MODELS
MU_HYDROPHOBIC = 4.0
SIGMA_HYDROPHOBIC = 1.5

MU_HYDROPHILIC = 3.0
SIGMA_HYDROPHILIC = 0.15

HPHOB_RNA_SUGAR = -0.13
HPHOB_RNA_PHOSPHATE = -2.02

MU_ANGLE_HBA = 129.9
MU_ANGLE_HBD = 109.0
MU_DIST_HBOND = 2.93

SIGMA_DIST_HBOND = 0.21
SIGMA_ANGLE_HBA = 20.0 # 29.3 originally (given precisely)
SIGMA_ANGLE_HBD = 20.0 # 20.0 originally (inferred from .figure)


MU_HBA = _np.array([MU_ANGLE_HBA, MU_DIST_HBOND])
COV_HBA = _np.array(
    [[SIGMA_ANGLE_HBA**2, 0],
    [0, SIGMA_DIST_HBOND**2]]
)
COV_INV_HBA = _np.linalg.inv(COV_HBA)


MU_HBD = _np.array([MU_ANGLE_HBD, MU_DIST_HBOND])
COV_HBD =  _np.array(
    [[SIGMA_ANGLE_HBD**2, 0],
    [0, SIGMA_DIST_HBOND**2]]
)
COV_INV_HBD = _np.linalg.inv(COV_HBD)


MU_STACKING = _np.array([29.97767535, 4.1876158])
COV_STACKING = _np.array(
    [[169.9862228 , 6.62318852],
     [  6.62318852, 0.37123882]]
)
COV_INV_STACKING = _np.linalg.inv(COV_STACKING)


### square root of the DIST contribution to COV_STACKING,
### calculated with _np.sqrt(COV_STACKING[1,1])
SIGMA_DIST_STACKING = 0.6092937058594976

######################## KERNEL_PARAMETERS

GAUSSIAN_KERNEL_SIGMAS = 4 # how many sigmas of width should the precalculated gaussians have?


################################################################################
################################### MODULES ####################################
################################################################################

from .core.grid import Grid, GridSMIF
from .core.kernel import Kernel
from .core.smiffer import Smiffer
from .core.systems import MolecularSystem, MSPocketSphere, MSWhole
from .core.trimmers import GridTrimmer, TrimmerPocketSphere, TrimmerWhole

from .grids.apbs import GridAPBS
# from .grids.cavities import GridCavities
from .grids.hbonds import GridHBAccepts, GridHBDonors
from .grids.hydro import GridHydrophilic, GridHydrophobic
from .grids.stacking import GridStacking

from .kernels.boolean import \
    KernelSphere, KernelCylinder, KernelDisk, KernelDiskConecut
from .kernels.gaussian import \
    KernelGaussianUnivariate, KernelGaussianMultivariate

from .utils.args import SmifferArgsParser
from .utils.io import \
    read_mrc, read_dx, read_cmap, \
    write_mrc, write_dx, write_cmap, \
    read_auto, grid_init_metadata, get_cmap_keys
from .utils.math import \
    normalize, dot_product, get_norm, get_angle, \
    get_projection, get_projection_height, \
    univariate_gaussian, multivariate_gaussian, \
    interpolate_3d, format_vector_str, get_coords_array
from .utils.tables import \
    planar_prot, planar_rna, ww_scale, \
    nucleic_backbone_phosphate, nucleic_backbone_sugar, \
    nucleic_bases, prot_hba, prot_hbd, rna_hba, rna_hbd
from .utils.timer import Timer


################################################################################
