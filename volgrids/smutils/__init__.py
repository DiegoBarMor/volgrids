from ._operations.rna_resids import RNAResids
from ._operations.dist import Histogram

from ._occupancies.hbaccepts import OgHBAccepts
from ._occupancies.hbdonors import OgHBDonors
from ._occupancies.stacking import OgStacking
from ._occupancies.hydrophobic import OgHydrophobic
# [TODO] hydrophilic? electrostatic?

from ._occupancies.app_occupancy import AppOccupancy
from .app_smutils import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

RADIUS_OCCUPANCY_OG: float = 2.0

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")
