from ._core.rna_resids import RNAResids
from ._core.operations import SMOperations

from ._core.app_smutils import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

### PLACEHOLDER: default configs go here

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")
