from pathlib import Path

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
        if operation == "convert"  : return self._run_convert()
        if operation == "pack"     : return self._run_pack()
        if operation == "unpack"   : return self._run_unpack()
        if operation == "fix_cmap" : return self._run_fix_cmap()
        if operation == "rotate"   : return self._run_rotate()
        if operation == "average"  : return self._run_average()
        if operation == "op"       : return self._run_op()
        if operation == "summary"  : return self._run_summary()
        if operation == "histogram": return self._run_histogram()
        if operation == "compare"  : return self._run_compare()
        if operation == "points"   : return self._run_points()


    # --------------------------------------------------------------------------
    def _run_convert(self):
        def _handle_out_arg(key: str, default_suffix: str) -> Path | None:
            path = self.main.get_arg_path(key, allow_none = True)

            #### CASE 0: flag not provided -> do not convert to this format
            if path is None: return

            #### CASE 1: flag provided with no path -> convert to this format, using input path with new suffix
            if path is True: return path_in.with_suffix(default_suffix)

            #### CASE 2: flag provided with path -> convert to this format, using provided path
            err = fy.PathAssertion.FILE_OUT(path)
            if isinstance(err, fy.ArgDTypeError):
                self.main.help_and_exit(1, err.err_message)
            return path

        def _convert(path_out, fmt_out: vg.GridFormat):
            if path_out is None: return
            print(f">>> Converting {fy.Color.yellow(path_in)} file to {fy.Color.magenta(fmt_out.name)}: {fy.Color.blue(path_out)}")
            vgt.VGOperations.convert(path_in, path_out, fmt_out)

        path_in       = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_out_dx   = _handle_out_arg("out_dx",   ".dx")
        path_out_mrc  = _handle_out_arg("out_mrc",  ".mrc")
        path_out_ccp4 = _handle_out_arg("out_ccp4", ".ccp4")
        path_out_cmap = _handle_out_arg("out_cmap", ".cmap")

        _convert(path_out_dx,   vg.GridFormat.DX)
        _convert(path_out_mrc,  vg.GridFormat.MRC)
        _convert(path_out_ccp4, vg.GridFormat.CCP4)
        _convert(path_out_cmap, vg.GridFormat.CMAP)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        paths_in = self.main.get_arg_path("paths_in", assertion = fy.PathAssertion.FILE_IN, is_list = True)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        str_npaths = fy.Color.yellow(f"{len(paths_in)}")
        print(f">>> Packing {str_npaths} grids into '{fy.Color.blue(path_out)}'")
        vgt.VGOperations.pack(paths_in, path_out)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("folder_out",
            assertion = fy.PathAssertion.DIR_OUT, default = path_in.parent
        )
        print(f">>> Unpacking '{fy.Color.yellow(path_in)}' into '{fy.Color.blue(path_out)}'")
        vgt.VGOperations.unpack(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_fix_cmap(self):
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        print(f">>> Fixing CMAP file: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.fix_cmap(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_rotate(self):
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        rotate_yz = self.main.get_arg_float("x")
        rotate_xz = self.main.get_arg_float("y")
        rotate_xy = self.main.get_arg_float("z")

        print(f">>> Rotating grid: {fy.Color.yellow(path_in)} by {rotate_xy}° (xy), {rotate_yz}° (yz), {rotate_xz}° (xz)")
        grid = vgt.VGOperations.rotate(path_in, rotate_xy, rotate_yz, rotate_xz)
        vg.GridIO.write_auto(path_out, grid)


    # --------------------------------------------------------------------------
    def _run_average(self):
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        print(f">>> Averaging CMAP file: {fy.Color.yellow(path_in)}")
        grid = vgt.VGOperations.average(path_in)
        vg.GridIO.write_auto(path_out, grid)


    # --------------------------------------------------------------------------
    def _run_op(self):
        command = self.main.subcommands.pop(0)

        if command == "abs": # abs is the only unary operation for now
            self._run_op_unary(command)
            return

        path_in_0 = self.main.get_arg_path("path_in_0", assertion = fy.PathAssertion.FILE_IN)
        path_in_1 = self.main.get_arg_path("path_in_1", assertion = fy.PathAssertion.FILE_IN)
        path_out  = self.main.get_arg_path("path_out",  assertion = fy.PathAssertion.FILE_OUT)

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
        for key, grid in vgt.VGOperations.iter_op_binary(
            path_in_0, path_in_1, operation, interpolate_to_common_box
        ):
            vg.GridIO.write_auto(path_out, grid, key)


    # --------------------------------------------------------------------------
    def _run_op_unary(self, command: str):
        operation: callable = {
            "abs": vg.Grid.__abs__,
        }[command]

        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        print(f">>> Performing '{fy.Color.yellow(command)}' operation on grid: {fy.Color.red(path_in)}")

        for key, grid in vgt.VGOperations.iter_op_unary(path_in, operation):
            vg.GridIO.write_auto(path_out, grid, key)


    # --------------------------------------------------------------------------
    def _run_summary(self):
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)

        print(f">>> Grid summary: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.summary(path_in)


    # --------------------------------------------------------------------------
    def _run_histogram(self) -> None:
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT, allow_none = True)
        key      = self.main.get_arg_str("key")

        print(f">>> Voxel distribution: {fy.Color.yellow(path_in)}")
        vgt.Histogram.plot(path_in, path_out, key)


    # --------------------------------------------------------------------------
    def _run_compare(self):
        path_in_0 = self.main.get_arg_path("path_0", assertion = fy.PathAssertion.FILE_IN)
        path_in_1 = self.main.get_arg_path("path_1", assertion = fy.PathAssertion.FILE_IN)

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
    def _run_points(self):
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        points_flat = self.main.get_arg_float("points", is_list = True)
        if len(points_flat) % 3: self.main.help_and_exit(1,
            f"Points should be provided as a list of floats, with 3 floats per XYZ point. "+\
            f"Got {len(points_flat)} floats, which is not a multiple of 3."
        )

        points = zip(*[points_flat[i::3] for i in range(3)])
        print(*vgt.VGOperations.points(path_in, *points))


# //////////////////////////////////////////////////////////////////////////////
