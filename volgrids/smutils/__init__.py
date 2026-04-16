from ._core.rna_resids import RNAResids
from ._core.operations import SMOperations

from ._ui.param_handler import ParamHandlerSMUtils
from ._ui.app import AppSMUtils


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

### PLACEHOLDER: default configs go here

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

OPERATION: str = '' # mode of the application, i.e. "resids_nonbp",

import pathlib as _pathlib

### (Generic)
PATH_STRUCT: _pathlib.Path = None # "path/input/structure.pdb"
