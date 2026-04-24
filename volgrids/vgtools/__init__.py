from ._operations.operations import VGOperations

from ._operations.comparison_result import ComparisonResult

from .app_vgtools import AppVGTools


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

### PLACEHOLDER: default configs go here

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")
