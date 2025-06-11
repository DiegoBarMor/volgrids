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
SAVE_CACHED_MASK = False # saves the logical inverse of the trimming mask

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

from .grids.apbs import GridAPBS
# from .grids.cavities import GridCavities
from .grids.hbonds import GridHBAccepts, GridHBDonors
from .grids.hydro import GridHydrophilic, GridHydrophobic
from .grids.stacking import GridStacking

from .args import SmifferArgsParser
from .tables import \
    planar_prot, planar_rna, ww_scale, \
    nucleic_backbone_phosphate, nucleic_backbone_sugar, \
    nucleic_bases, prot_hba, prot_hbd, rna_hba, rna_hbd
from .trimmers import GridTrimmer


from .calculator import SmifferCalculator
