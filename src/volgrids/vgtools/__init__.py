from ._core.operations import VGOperations

from ._misc.comparison_result import ComparisonResult

from ._ui.param_handler import ParamHandlerVGTools
from ._ui.app import AppVGTools


__config_keys__ = set(__annotations__.keys())


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

OPERATION: str = '' # mode of the application, i.e. "convert", "pack", "unpack", "fix_cmap"

import pathlib as _pathlib

### Convert
PATH_CONVERT_IN:   _pathlib.Path = None # "path/input/grid.dx" or "path/input/grid.mrc" or "path/input/grid.cmap"
PATH_CONVERT_DX:   _pathlib.Path = None # "path/input/grid.dx"
PATH_CONVERT_MRC:  _pathlib.Path = None # "path/output/grid.mrc"
PATH_CONVERT_CCP4: _pathlib.Path = None # "path/output/grid.ccp4"
PATH_CONVERT_CMAP: _pathlib.Path = None # "path/output/grid.cmap"

### Pack
PATHS_PACK_IN: list[_pathlib.Path] = None # list of paths to input grids for packing
PATH_PACK_OUT: _pathlib.Path = None # "path/output/packed.cmap"

### Unpack
PATH_UNPACK_IN:  _pathlib.Path = None # "path/input/packed.cmap"
PATH_UNPACK_OUT: _pathlib.Path = None # folder where to unpack the grids

### Fix CMAP
PATH_FIXCMAP_IN:  _pathlib.Path = None # "path/input/fix.cmap"
PATH_FIXCMAP_OUT: _pathlib.Path = None # "path/output/fix.cmap"

### Compare
PATH_COMPARE_IN_0: _pathlib.Path = None # "path/input/grid_0.mrc"
PATH_COMPARE_IN_1: _pathlib.Path = None # "path/input/grid_1.mrc"
THRESHOLD_COMPARE: float # threshold for comparison (default 1e-5)
