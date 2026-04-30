import volgrids as vg
import volgrids.vgtools as vgt
from volgrids._vendors import freyacli as fy

DEFAULT_COMPARISON_THRESHOLD = 1e-5 # [TODO] this should be a config

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.AppSubcommand):
    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        app_main.load_configs(vg, vgt)


    # --------------------------------------------------------------------------
    def run(self) -> None:
        operation = self.main.subcommands.pop(0)
        if operation == "convert" : return self._run_convert()
        if operation == "pack"    : return self._run_pack()
        if operation == "unpack"  : return self._run_unpack()
        if operation == "fix_cmap": return self._run_fix_cmap()
        if operation == "average" : return self._run_average()
        if operation == "summary" : return self._run_summary()
        if operation == "compare" : return self._run_compare()
        if operation == "rotate"  : return self._run_rotate()
        if operation == "op"      : return self._run_op()


    # --------------------------------------------------------------------------
    def _run_convert(self):
        path_in       = self.main.get_arg_path("path_in")
        path_out_dx   = self.main.get_arg_path("out_dx")
        path_out_mrc  = self.main.get_arg_path("out_mrc")
        path_out_ccp4 = self.main.get_arg_path("out_ccp4")
        path_out_cmap = self.main.get_arg_path("out_cmap")

        if path_out_dx is True: # -d flag with no path specified
            path_out_dx = path_in.with_suffix(".dx")

        if path_out_mrc is True: # -m flag with no path specified
            path_out_mrc = path_in.with_suffix(".mrc")

        if path_out_ccp4 is True: # -p flag with no path specified
            path_out_ccp4 = path_in.with_suffix(".ccp4")

        if path_out_cmap is True: # -c flag with no path specified
            path_out_cmap = path_in.with_suffix(".cmap")

        self.main.assert_file_in(path_in)
        self.main.assert_file_out(path_out_dx, allow_none = True)
        self.main.assert_file_out(path_out_mrc, allow_none = True)
        self.main.assert_file_out(path_out_ccp4, allow_none = True)
        self.main.assert_file_out(path_out_cmap, allow_none = True)


        def _convert(path_out, fmt_out: vg.GridFormat):
            if path_out is None: return
            print(f">>> Converting {fy.Color.yellow(path_in)} file to {fy.Color.magenta(fmt_out.name)}: {fy.Color.blue(path_out)}")
            vgt.VGOperations.convert(path_in, path_out, fmt_out)

        _convert(path_out_dx,   vg.GridFormat.DX)
        _convert(path_out_mrc,  vg.GridFormat.MRC)
        _convert(path_out_ccp4, vg.GridFormat.CCP4)
        _convert(path_out_cmap, vg.GridFormat.CMAP)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        paths_in = self.main.get_arg_list_path("paths_in")
        path_out = self.main.get_arg_path("path_out")

        for path in paths_in: self.main.assert_file_in(path)
        self.main.assert_file_out(path_out)

        str_npaths = fy.Color.yellow(f"{len(paths_in)}")
        print(f">>> Packing {str_npaths} grids into '{fy.Color.blue(path_out)}'")
        vgt.VGOperations.pack(paths_in, path_out)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("folder_out", default = path_in.parent)

        self.main.assert_file_in(path_in)
        self.main.assert_dir_out(path_out)

        print(f">>> Unpacking '{fy.Color.yellow(path_in)}' into '{fy.Color.blue(path_out)}'")
        vgt.VGOperations.unpack(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_fix_cmap(self):
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")

        self.main.assert_file_in(path_in)
        self.main.assert_file_out(path_out)

        print(f">>> Fixing CMAP file: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.fix_cmap(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_average(self):
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")

        self.main.assert_file_in(path_in)
        self.main.assert_dir_out(path_out)

        print(f">>> Averaging CMAP file: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.average(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_summary(self):
        path_in = self.main.get_arg_path("path_in")
        self.main.assert_file_in(path_in)

        print(f">>> Grid summary: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.summary(path_in)


    # --------------------------------------------------------------------------
    def _run_compare(self):
        path_in_0 = self.main.get_arg_path("path_0")
        path_in_1 = self.main.get_arg_path("path_1")

        self.main.assert_file_in(path_in_0)
        self.main.assert_file_in(path_in_1)

        threshold = self.main.get_arg_float("threshold", default = DEFAULT_COMPARISON_THRESHOLD)

        print(f">>> Comparing grids: {fy.Color.red(path_in_0)} vs {fy.Color.blue(path_in_1)} (threshold={threshold:2.2e})")
        result = vgt.VGOperations.compare(path_in_0, path_in_1, threshold)

        for message in result.messages:
            print(f"...>>> {message}")
        if result.npoints_total == 0: return

        str_perc = fy.Color.yellow(f"{100 * result.npoints_diff / result.npoints_total:.2f}%")
        str_avg_diff = fy.Color.red(f"{result.avg_diff:2.2e}")

        print(
            f"...>>> {result.npoints_diff}/{result.npoints_total} points differ ({str_perc})\n" +\
            f"...>>> Accumulated difference: {result.cumulative_diff:2.2e} (avg {str_avg_diff} per point)"
        )


    # --------------------------------------------------------------------------
    def _run_rotate(self):
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")

        self.main.assert_file_in(path_in)
        self.main.assert_dir_out(path_out)

        rotate_yz = self.main.get_arg_float("x")
        rotate_xz = self.main.get_arg_float("y")
        rotate_xy = self.main.get_arg_float("z")

        print(f">>> Rotating grid: {fy.Color.yellow(path_in)} by {rotate_xy}° (xy), {rotate_yz}° (yz), {rotate_xz}° (xz)")
        vgt.VGOperations.rotate(path_in, path_out, rotate_xy, rotate_yz, rotate_xz)


    # --------------------------------------------------------------------------
    def _run_op(self):
        command = self.main.subcommands.pop(0)

        path_out  = self.main.get_arg_path("path_out")
        self.main.assert_file_out(path_out)

        if command == "abs": # abs is the only unary operation for now
            operation = vg.Grid.__abs__
            path_in = self.main.get_arg_path("path_in")
            self.main.assert_file_in(path_in)
            print(f">>> Performing '{fy.Color.yellow(command)}' operation on grid: {fy.Color.red(path_in)}")
            vgt.VGOperations.op(operation, path_out, path_in)
            return

        path_in_0 = self.main.get_arg_path("path_in_0")
        path_in_1 = self.main.get_arg_path("path_in_1")

        self.main.assert_file_in(path_in_0)
        self.main.assert_file_in(path_in_1)

        interpolate_to_common_box = self.main.get_arg_bool("common_box")

        operation: callable = {
            "add": vg.Grid.__add__,
            "sub": vg.Grid.__sub__,
            "mul": vg.Grid.__mul__,
            "div": vg.Grid.__div__,
            "and": vg.Grid.__and__,
            "or" : vg.Grid.__or__,
        }[command]

        print(f">>> Performing '{fy.Color.yellow(command)}' operation on grids: {fy.Color.red(path_in_0)} vs {fy.Color.blue(path_in_1)}")
        vgt.VGOperations.op(operation, path_out, path_in_0, path_in_1, interpolate_to_common_box)


# //////////////////////////////////////////////////////////////////////////////
