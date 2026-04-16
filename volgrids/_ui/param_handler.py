import volgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerVGTools(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
            "help"      : ("-h", "--help"),
            "input"     : ("-i", "--input"),
            "output"    : ("-o", "--output"),
            "dx"        : ("-d", "--dx"),
            "mrc"       : ("-m", "--mrc"),
            "ccp4"      : ("-p", "--ccp4"),
            "cmap"      : ("-c", "--cmap"),
            "thresh"    : ("-t", "--threshold"),
            "rot_x"     : ("-x", "--yz"),
            "rot_y"     : ("-y", "--xz"),
            "rot_z"     : ("-z", "--xy"),
            "operation" : ("--op", "--operation"),
    }
    _DEFAULT_COMPARISON_THRESHOLD = 1e-5


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: volgrids vgtools [convert|pack|unpack|average|fix_cmap|overlap|overlap_cross|overlap_diff] [options...]",
            "Available modes:",
            "    convert      - Convert grid files between formats.",
            "    pack         - Pack multiple grid files into a single CMAP series-file.",
            "    unpack       - Unpack a CMAP series-file into multiple grid files.",
            "    average      - Average all grids in a CMAP series-file into a single grid.",
            "    fix_cmap     - Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.",
            "    summary      - Print a summary of the grid file (format, dimensions, resolution, etc.) to the console.",
            "    compare      - Compare two grid files by printing the number of differing points and their accumulated difference.",
            "    rotate       - Rotate a grid file by 3 angles, along the xy, yz and xz planes (in degrees).",
            "    overlap      - Compute overlap between two molecular interaction fields.",
            "    overlap_cross - Smart cross-comparison overlap analysis between different field types.",
            "    overlap_diff  - Compute difference grids between two molecular interaction fields.",
            "\nRun 'volgrids vgtools [mode] --help' for more details on each mode.",
        )
        if self._has_param_kwds("help") and not self._has_params_pos():
            self._exit_with_help()

        vgt.OPERATION = self._safe_get_param_pos(0).lower()
        func: callable = self._safe_map_value(vgt.OPERATION,
            convert       = self._parse_convert,
            pack          = self._parse_pack,
            unpack        = self._parse_unpack,
            fix_cmap      = self._parse_fix_cmap,
            average       = self._parse_average,
            summary       = self._parse_summary,
            compare       = self._parse_compare,
            rotate        = self._parse_rotate,
            overlap       = self._parse_overlap,
            overlap_cross = self._parse_overlap_cross,
            overlap_diff  = self._parse_overlap_diff,
        )
        func()


    # --------------------------------------------------------------------------
    def _parse_convert(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools convert [input_grid] [options...]",
            "Available options:",
            "    -h, --help  Show this help message and exit.",
            "    -d, --dx    File path where to save the converted grid in DX format.",
            "    -m, --mrc   File path where to save the converted grid in MRC format.",
            "    -p, --ccp4  File path where to save the converted grid in CCP4 format.",
            "    -c, --cmap  File path where to save the converted grid in CMAP format. The stem of the input file will be used as the CMAP key.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_CONVERT_IN = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input grid file provided. Provide a path to the grid file as second positional argument."
            )
        )

        vgt.PATH_CONVERT_DX   = self._safe_kwd_file_out("dx")
        vgt.PATH_CONVERT_MRC  = self._safe_kwd_file_out("mrc")
        vgt.PATH_CONVERT_CCP4 = self._safe_kwd_file_out("ccp4")
        vgt.PATH_CONVERT_CMAP = self._safe_kwd_file_out("cmap")


    # --------------------------------------------------------------------------
    def _parse_pack(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools pack [options...]",
            "Available options:",
            "    -h, --help    Show this help message and exit.",
            "    -i, --input   List of file paths with the input grids to be packed. At least one grid file must be provided.",
            "    -o, --output  File path where to save the packed grid in CMAP format. Must be provided.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATHS_PACK_IN = [
            self._safe_path_file_in(path) for path in \
            self._safe_get_param_kwd_list("input")
        ]

        vgt.PATH_PACK_OUT = self._safe_kwd_file_out("output", required = True)


    # --------------------------------------------------------------------------
    def _parse_unpack(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools unpack [options...]",
            "Available options:",
            "    -h, --help    Show this help message and exit.",
            "    -i, --input   File path to the CMAP series-file to be unpacked. Must be provided.",
            "    -o, --output  Folder path where to save the unpacked grids. If not provided, the parent folder of the input packed file will be used.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_UNPACK_IN  = self._safe_kwd_file_in("input", required = True)
        vgt.PATH_UNPACK_OUT = self._safe_kwd_folder_out("output", default = vgt.PATH_UNPACK_IN.parent)


    # --------------------------------------------------------------------------
    def _parse_fix_cmap(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools fix_cmap [options...]",
            "Available options:",
            "    -h, --help    Show this help message and exit.",
            "    -i, --input   File path to the CMAP file to be fixed. Must be provided.",
            "    -o, --output  File path where to save the fixed CMAP file. Must be provided.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_FIXCMAP_IN  = self._safe_kwd_file_in("input", required = True)
        vgt.PATH_FIXCMAP_OUT = self._safe_kwd_file_out("output", required = True)


    # --------------------------------------------------------------------------
    def _parse_average(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools average [options...]",
            "Available options:",
            "    -h, --help    Show this help message and exit.",
            "    -i, --input   File path to the CMAP series-file to be averaged. Must be provided.",
            "    -o, --output  File path where to save the averaged CMAP file. Must be provided.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_AVERAGE_IN  = self._safe_kwd_file_in("input", required = True)
        vgt.PATH_AVERAGE_OUT = self._safe_kwd_file_out("output", required = True)


    # --------------------------------------------------------------------------
    def _parse_summary(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools summary [input_grid] [options...]",
            "Available options:",
            "    -h, --help  Show this help message and exit.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_SUMMARY_IN = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input grid file provided. Provide a path to the grid file as second positional argument."
            )
        )


    # --------------------------------------------------------------------------
    def _parse_compare(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools compare [input_grid_0] [input_grid_1] [options...]",
            "Available options:",
            "    -h, --help       Show this help message and exit.",
            "    -t, --threshold  Threshold for comparison. Default is 1e-3.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_COMPARE_IN_0 = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No first input grid file provided. Provide a path to the first grid file as second positional argument."
            )
        )

        vgt.PATH_COMPARE_IN_1 = self._safe_path_file_in(
            self._safe_get_param_pos(2,
               err_msg = "No second input grid file provided. Provide a path to the second grid file as third positional argument."
            )
        )

        vgt.THRESHOLD_COMPARE = self._safe_kwd_float("thresh", default = self._DEFAULT_COMPARISON_THRESHOLD)


    # --------------------------------------------------------------------------
    def _parse_rotate(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools rotate [input_grid] [output_grid]",
            "Available options:",
            "    -h, --help  Show this help message and exit.",
            "    -x, --yz    Rotation angle in degrees along the yz plane (in degrees).",
            "    -y, --xz    Rotation angle in degrees along the xz plane (in degrees).",
            "    -z, --xy    Rotation angle in degrees along the xy plane (in degrees).",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_ROTATE_IN = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input grid file provided. Provide a path to the grid file as second positional argument."
            )
        )
        vgt.PATH_ROTATE_OUT = self._safe_path_file_out(
            self._safe_get_param_pos(2,
               err_msg = "No output grid file provided. Provide a path where to save the rotated grid as third positional argument."
            )
        )
        vgt.ROTATE_YZ = self._safe_kwd_float("rot_x", default = 0.0)
        vgt.ROTATE_XZ = self._safe_kwd_float("rot_y", default = 0.0)
        vgt.ROTATE_XY = self._safe_kwd_float("rot_z", default = 0.0)


    # --------------------------------------------------------------------------
    def _parse_overlap(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools overlap [grid1] [grid2] [output] [options...]",
            "Compute overlap between two molecular interaction fields.",
            "The first grid will be interpolated to match the coordinate system of the second grid.",
            "Available options:",
            "    -h, --help         Show this help message and exit.",
            "    --op, --operation  Operation type: 'multiply' (default), 'subtract', 'add'.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_OVERLAP_GRID1 = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No first grid file provided. Provide a path to the first grid file as second positional argument."
            )
        )

        vgt.PATH_OVERLAP_GRID2 = self._safe_path_file_in(
            self._safe_get_param_pos(2,
               err_msg = "No second grid file provided. Provide a path to the second grid file as third positional argument."
            )
        )

        vgt.PATH_OVERLAP_OUT = self._safe_path_file_out(
            self._safe_get_param_pos(3,
               err_msg = "No output grid file provided. Provide a path where to save the overlap grid as fourth positional argument."
            )
        )

        vgt.OVERLAP_OPERATION = self._safe_kwd_str("operation", default="multiply")


    # --------------------------------------------------------------------------
    def _parse_overlap_cross(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools overlap_cross [grid1] [grid2] [output]",
            "Cross-comparison overlap analysis (multiplication).",
            "Available options:",
            "    -h, --help  Show this help message and exit.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_OVERLAP_CROSS_GRID1 = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No first grid file provided. Provide a path to the first grid file as second positional argument."
            )
        )

        vgt.PATH_OVERLAP_CROSS_GRID2 = self._safe_path_file_in(
            self._safe_get_param_pos(2,
               err_msg = "No second grid file provided. Provide a path to the second grid file as third positional argument."
            )
        )

        vgt.PATH_OVERLAP_CROSS_OUT = self._safe_path_file_out(
            self._safe_get_param_pos(3,
               err_msg = "No output grid file provided. Provide a path where to save the overlap grid as fourth positional argument."
            )
        )


    # --------------------------------------------------------------------------
    def _parse_overlap_diff(self) -> None:
        self._set_help_str(
            "usage: volgrids vgtools overlap_diff [grid1] [grid2] [output]",
            "Difference overlap analysis (subtraction: grid1 - grid2).",
            "Available options:",
            "    -h, --help  Show this help message and exit.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        vgt.PATH_OVERLAP_DIFF_GRID1 = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No first grid file provided. Provide a path to the first grid file as second positional argument."
            )
        )

        vgt.PATH_OVERLAP_DIFF_GRID2 = self._safe_path_file_in(
            self._safe_get_param_pos(2,
               err_msg = "No second grid file provided. Provide a path to the second grid file as third positional argument."
            )
        )

        vgt.PATH_OVERLAP_DIFF_OUT = self._safe_path_file_out(
            self._safe_get_param_pos(3,
               err_msg = "No output grid file provided. Provide a path where to save the difference grid as fourth positional argument."
            )
        )


# //////////////////////////////////////////////////////////////////////////////
