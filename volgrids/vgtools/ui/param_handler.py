import volgrids.vgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerVGTools(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
            "help"   : ("-h", "--help"),
            "input"  : ("-i", "--input"),
            "output" : ("-o", "--output"),
            "dx"     : ("-d", "--dx"),
            "mrc"    : ("-m", "--mrc"),
            "ccp4"   : ("-p", "--ccp4"),
            "cmap"   : ("-c", "--cmap"),
    }


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: python3 vgtools.py [convert|pack|unpack|fix_cmap] [options...]",
            "Available modes:",
            "  convert  - Convert grid files between formats.",
            "  pack     - Pack multiple grid files into a single CMAP series-file.",
            "  unpack   - Unpack a CMAP series-file into multiple grid files.",
            "  fix_cmap - Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.",
            "Run 'python3 vgtools.py [mode] --help' for more details on each mode.",
        )
        if self._has_param_kwds("help") and not self._has_params_pos():
            self._exit_with_help(0)

        vgt.OPERATION = self._safe_get_param_pos(0).lower()
        func: callable = self._safe_map_value(vgt.OPERATION,
            convert  = self._parse_convert,
            pack     = self._parse_pack,
            unpack   = self._parse_unpack,
            fix_cmap = self._parse_fix_cmap,
        )
        func()


    # --------------------------------------------------------------------------
    def _parse_convert(self) -> None:
        self._set_help_str(
            "usage: python3 vgtools.py convert [path/input/grid.dx] [options...]",
            "Available options:",
            "-h, --help  Show this help message and exit.",
            "-d, --dx    File path where to save the converted grid in DX format.",
            "-m, --mrc   File path where to save the converted grid in MRC format.",
            "-p, --ccp4  File path where to save the converted grid in CCP4 format.",
            "-c, --cmap  File path where to save the converted grid in CMAP format. The stem of the input file will be used as the CMAP key.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help(0)

        vgt.PATH_CONVERT_IN = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input grid file provided. Provide a path to the grid file as first positional argument."
            )
        )

        if self._has_param_kwds("dx"):
            vgt.PATH_CONVERT_DX = self._safe_kwd_file_out("dx")

        if self._has_param_kwds("mrc"):
            vgt.PATH_CONVERT_MRC = self._safe_kwd_file_out("mrc")

        if self._has_param_kwds("ccp4"):
            vgt.PATH_CONVERT_CCP4 = self._safe_kwd_file_out("ccp4")

        if self._has_param_kwds("cmap"):
            vgt.PATH_CONVERT_CMAP = self._safe_kwd_file_out("cmap")


    # --------------------------------------------------------------------------
    def _parse_pack(self) -> None:
        self._set_help_str(
            "usage: python3 vgtools.py pack [options...]",
            "Available options:",
            "-h, --help    Show this help message and exit.",
            "-i, --input   List of file paths with the input grids to be packed. At least one grid file must be provided.",
            "-o, --output  File path where to save the packed grid in CMAP format. Must be provided.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help(0)

        vgt.PATHS_PACK_IN = [
            self._safe_path_file_in(path) for path in \
            self._safe_get_param_kwd_list("input", 0)
        ]

        vgt.PATH_PACK_OUT = self._safe_kwd_file_out("output")


    # --------------------------------------------------------------------------
    def _parse_unpack(self) -> None:
        self._set_help_str(
            "usage: python3 vgtools.py unpack [options...]",
            "Available options:",
            "-h, --help    Show this help message and exit.",
            "-i, --input   File path to the CMAP series-file to be unpacked. Must be provided.",
            "-o, --output  Folder path where to save the unpacked grids. If not provided, the parent folder of the input packed file will be used.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help(0)

        vgt.PATH_UNPACK_IN = self._safe_kwd_file_in("input")

        vgt.PATH_UNPACK_OUT = self._safe_path_folder_out(
            self._safe_get_param_kwd("output", 0) if self._has_param_kwds("output") \
            else vgt.PATH_UNPACK_IN.parent
        )


    # --------------------------------------------------------------------------
    def _parse_fix_cmap(self) -> None:
        self._set_help_str(
            "usage: python3 vgtools.py fix_cmap [options...]",
            "Available options:",
            "-h, --help    Show this help message and exit.",
            "-i, --input   File path to the CMAP file to be fixed. Must be provided.",
            "-o, --output  File path where to save the fixed CMAP file. Must be provided.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help(0)

        vgt.PATH_FIXCMAP_IN = self._safe_kwd_file_in("input")

        vgt.PATH_FIXCMAP_OUT = self._safe_kwd_file_out("output")


# //////////////////////////////////////////////////////////////////////////////
