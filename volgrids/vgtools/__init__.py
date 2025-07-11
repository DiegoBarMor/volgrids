from .ui.param_handler import ParamHandlerVGTools

from .vgtools import VGToolsApp


######################## COMMAND LINE ARGUMENTS GLOBALS ########################
### These are global variables that are to be set by
### an instance of ParamHandler (or its inherited classes)

from pathlib import Path

### Convert
PATH_CONVERT_IN:   Path = None # "path/input/grid.dx" or "path/input/grid.mrc" or "path/input/grid.cmap"
PATH_CONVERT_DX:   Path = None # "path/input/grid.dx"
PATH_CONVERT_MRC:  Path = None # "path/output/grid.mrc"
PATH_CONVERT_CCP4: Path = None # "path/output/grid.ccp4"
PATH_CONVERT_CMAP: Path = None # "path/output/grid.cmap"

### Pack
PATHS_PACK_IN: list[Path] = None # list of paths to input grids for packing
PATH_PACK_OUT: Path = None # "path/output/packed.cmap"

### Unpack
PATH_UNPACK_IN:  Path = None # "path/input/packed.cmap"
PATH_UNPACK_OUT: Path = None # folder where to unpack the grids

### Fix CMAP
PATH_FIXCMAP_IN:  Path = None # "path/input/fix.cmap"
PATH_FIXCMAP_OUT: Path = None # "path/output/fix.cmap"
