import volgrids as vg
import volgrids.vgtools as vgt

### [TODO] many of the globals here could be refactored into local variables
DEFAULT_COMPARISON_THRESHOLD = 1e-5 # [TODO] this should be a config

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.AppSubcommand):
    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.func_operation: callable = None

    # --------------------------------------------------------------------------
    def assign_globals(self):
        vgt.OPERATION = self.main.subcommands.pop(0)
        self.func_operation = {
            "convert" : self._run_convert,
            "pack"    : self._run_pack,
            "unpack"  : self._run_unpack,
            "fix_cmap": self._run_fix_cmap,
            "average" : self._run_average,
            "summary" : self._run_summary,
            "compare" : self._run_compare,
            "rotate"  : self._run_rotate,
        }[vgt.OPERATION]


    # --------------------------------------------------------------------------
    def run(self) -> None:
        self.func_operation()


    # --------------------------------------------------------------------------
    def _run_convert(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["out_dx", "out_mrc", "out_ccp4", "out_cmap"],
            allow_none = True,
        )
        vgt.PATH_CONVERT_IN   = self.main.get_arg_path("path_in")
        vgt.PATH_CONVERT_DX   = self.main.get_arg_path("out_dx")
        vgt.PATH_CONVERT_MRC  = self.main.get_arg_path("out_mrc")
        vgt.PATH_CONVERT_CCP4 = self.main.get_arg_path("out_ccp4")
        vgt.PATH_CONVERT_CMAP = self.main.get_arg_path("out_cmap")

        def _convert(path_out, fmt_out: vg.GridFormat):
            if path_out is None: return
            print(f">>> Converting {vgt.PATH_CONVERT_IN} file to {fmt_out.name}: {path_out}")
            vgt.VGOperations.convert(vgt.PATH_CONVERT_IN, path_out, fmt_out)

        _convert(vgt.PATH_CONVERT_DX,   vg.GridFormat.DX)
        _convert(vgt.PATH_CONVERT_MRC,  vg.GridFormat.MRC)
        _convert(vgt.PATH_CONVERT_CCP4, vg.GridFormat.CCP4)
        _convert(vgt.PATH_CONVERT_CMAP, vg.GridFormat.CMAP)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        self.main.assert_paths(
            keys_file_in = ["paths_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        vgt.PATHS_PACK_IN = self.main.get_arg_list_path("paths_in")
        vgt.PATH_PACK_OUT = self.main.get_arg_path("path_out")

        print(f">>> Packing {len(vgt.PATHS_PACK_IN)} grids into '{vgt.PATH_PACK_OUT}'")
        vgt.VGOperations.pack(vgt.PATHS_PACK_IN, vgt.PATH_PACK_OUT)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        self.main.assert_paths(keys_file_in = ["path_in"], allow_none = False)
        vgt.PATH_UNPACK_IN = self.main.get_arg_path("path_in")

        if self.main.get_arg_value("folder_out") is None:
            self.main.set_arg_value("folder_out", vgt.PATH_UNPACK_IN.parent)

        self.main.assert_paths(keys_dir_out = ["folder_out"], allow_none = True)
        vgt.PATH_UNPACK_OUT = self.main.get_arg_path("folder_out")

        print(f">>> Unpacking '{vgt.PATH_UNPACK_IN}' into '{vgt.PATH_UNPACK_OUT}'")
        vgt.VGOperations.unpack(vgt.PATH_UNPACK_IN, vgt.PATH_UNPACK_OUT)


    # --------------------------------------------------------------------------
    def _run_fix_cmap(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        vgt.PATH_FIXCMAP_IN  = self.main.get_arg_path("path_in")
        vgt.PATH_FIXCMAP_OUT = self.main.get_arg_path("path_out")

        print(f">>> Fixing CMAP file: {vgt.PATH_FIXCMAP_IN}")
        vgt.VGOperations.fix_cmap(vgt.PATH_FIXCMAP_IN, vgt.PATH_FIXCMAP_OUT)


    # --------------------------------------------------------------------------
    def _run_average(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        vgt.PATH_AVERAGE_IN  = self.main.get_arg_path("path_in")
        vgt.PATH_AVERAGE_OUT = self.main.get_arg_path("path_out")

        print(f">>> Averaging CMAP file: {vgt.PATH_AVERAGE_IN}")
        vgt.VGOperations.average(vgt.PATH_AVERAGE_IN, vgt.PATH_AVERAGE_OUT)


    # --------------------------------------------------------------------------
    def _run_summary(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            allow_none = False,
        )
        vgt.PATH_SUMMARY_IN = self.main.get_arg_path("path_in")

        print(f">>> Grid summary: {vgt.PATH_SUMMARY_IN}")
        vgt.VGOperations.summary(vgt.PATH_SUMMARY_IN)


    # --------------------------------------------------------------------------
    def _run_compare(self):
        self.main.assert_paths(
            keys_file_in = ["path_0", "path_1"],
            allow_none = False,
        )
        vgt.PATH_COMPARE_IN_0 = self.main.get_arg_path("path_0")
        vgt.PATH_COMPARE_IN_1 = self.main.get_arg_path("path_1")

        vgt.THRESHOLD_COMPARE = self.main.get_arg_float("threshold")
        if vgt.THRESHOLD_COMPARE is None:
            vgt.THRESHOLD_COMPARE = DEFAULT_COMPARISON_THRESHOLD

        print(f">>> Comparing grids: {vgt.PATH_COMPARE_IN_0} vs {vgt.PATH_COMPARE_IN_1} (threshold={vgt.THRESHOLD_COMPARE:2.2e})")
        result = vgt.VGOperations.compare(vgt.PATH_COMPARE_IN_0, vgt.PATH_COMPARE_IN_1, vgt.THRESHOLD_COMPARE)

        for message in result.messages:
            print(f"...>>> {message}")
        if result.npoints_total == 0: return

        print(
            f"...>>> {result.npoints_diff}/{result.npoints_total} points differ " +\
            f"({100 * result.npoints_diff / result.npoints_total:.2f}%)\n" +\
            f"...>>> Accumulated difference: {result.cumulative_diff:2.2e} " +\
            f"(avg {result.avg_diff:2.2e} per point)"
        )


    # --------------------------------------------------------------------------
    def _run_rotate(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        vgt.PATH_ROTATE_IN  = self.main.get_arg_path("path_in")
        vgt.PATH_ROTATE_OUT = self.main.get_arg_path("path_out")

        vgt.ROTATE_YZ = self.main.get_arg_float("x")
        vgt.ROTATE_XZ = self.main.get_arg_float("y")
        vgt.ROTATE_XY = self.main.get_arg_float("z")

        print(f">>> Rotating grid: {vgt.PATH_ROTATE_IN} by {vgt.ROTATE_XY}° (xy), {vgt.ROTATE_YZ}° (yz), {vgt.ROTATE_XZ}° (xz)")
        vgt.VGOperations.rotate(vgt.PATH_ROTATE_IN, vgt.PATH_ROTATE_OUT, vgt.ROTATE_XY, vgt.ROTATE_YZ, vgt.ROTATE_XZ)


# # //////////////////////////////////////////////////////////////////////////////
