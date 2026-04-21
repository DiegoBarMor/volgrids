from ._core.operations import VGOperations

from ._misc.comparison_result import ComparisonResult

from ._core.app_vgtools import AppVGTools


############################# CONFIG FILE GLOBALS ##############################
_keys_other = set(globals().keys())

### PLACEHOLDER: default configs go here

__config_keys__ = set(globals().keys()) - _keys_other
__config_keys__.remove("_keys_other")

### Overlap
PATH_OVERLAP_GRID1: _pathlib.Path = None # "path/input/grid1.dx" (first grid, will be interpolated)
PATH_OVERLAP_GRID2: _pathlib.Path = None # "path/input/grid2.dx" (second grid, defines target coordinate system)
PATH_OVERLAP_OUT:   _pathlib.Path = None # "path/output/overlap.dx" (output overlap grid)
OVERLAP_OPERATION:  str = "multiply"     # operation type: "multiply", "subtract", "add"

### Overlap Cross-Comparison
PATH_OVERLAP_CROSS_GRID1: _pathlib.Path = None # "path/input/grid1.dx"
PATH_OVERLAP_CROSS_GRID2: _pathlib.Path = None # "path/input/grid2.dx"
PATH_OVERLAP_CROSS_OUT:   _pathlib.Path = None # "path/output/cross_overlap.dx"

### Overlap Difference
PATH_OVERLAP_DIFF_GRID1: _pathlib.Path = None # "path/input/grid1.dx" (minuend)
PATH_OVERLAP_DIFF_GRID2: _pathlib.Path = None # "path/input/grid2.dx" (subtrahend)
PATH_OVERLAP_DIFF_OUT:   _pathlib.Path = None # "path/output/diff_overlap.dx"
