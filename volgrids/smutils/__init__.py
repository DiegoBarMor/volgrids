from ._core.rna_resids import RNAResids
from ._core.operations import SMOperations

from ._occupancies.hbaccepts import OgHBAccepts
from ._occupancies.hbdonors import OgHBDonors
from ._occupancies.stacking import OgStacking
from ._occupancies.hydrophobic import OgHydrophobic
# [TODO] hydrophilic? electrostatic?

from ._occupancies.app_occupancy import AppOccupancy
from ._core.app_smutils import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

### PLACEHOLDER: default configs go here

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")
