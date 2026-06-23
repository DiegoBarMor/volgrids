from pathlib import Path

import volgrids as vg
import volgrids.vgtools as vgt
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        app_main.load_configs(vg, vgt)


    # --------------------------------------------------------------------------
    def run(self) -> None:
        operation = self.main.subcommands.pop(0)
        if operation == "convert"  : return self._run_convert()
        if operation == "pack"     : return self._run_pack()
        if operation == "unpack"   : return self._run_unpack()
        if operation == "extract"  : return self._run_extract()
        if operation == "fix_cmap" : return self._run_fix_cmap()
        if operation == "rotate"   : return self._run_rotate()
        if operation == "segment"  : return self._run_segment()
        if operation == "average"  : return self._run_average()
        if operation == "std_dev"  : return self._run_std_dev()
        if operation == "op"       : return self._run_op()
        if operation == "summary"  : return self._run_summary()
        if operation == "histogram": return self._run_histogram()
        if operation == "compare"  : return self._run_compare()
        if operation == "points"   : return self._run_points()
        raise ValueError(f"Unknown operation: {operation}")


    # --------------------------------------------------------------------------
    def _run_convert(self):
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        folder_out = self.main.get_arg_path("folder_out",
            assertion = fy.PathAssertion.DIR_OUT, default = path_in.parent
        )
        fmt = self._get_valid_fmt_arg()

        path_out = folder_out / f"{path_in.stem}.{fmt.suffix()}"

        print(f">>> Converting {fy.Color.yellow(path_in)} file to {fy.Color.magenta(fmt.name)}: {fy.Color.blue(path_out)}")
        vgt.VGOperations.convert(path_in, path_out, fmt)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        paths_in = self.main.get_arg_path("paths_in", assertion = fy.PathAssertion.FILE_IN, is_list = True)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        str_npaths = fy.Color.yellow(f"{len(paths_in)}")
        print(f">>> Packing {str_npaths} grids into '{fy.Color.blue(path_out)}'")
        vgt.VGOperations.pack(paths_in, path_out)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        folder_out = self.main.get_arg_path("folder_out",
            assertion = fy.PathAssertion.DIR_OUT, default = path_in.parent
        )
        fmt = self._get_valid_fmt_arg()

        print(f">>> Unpacking '{fy.Color.yellow(path_in)}' into '{fy.Color.blue(folder_out)}'")
        vgt.VGOperations.unpack(path_in, folder_out, fmt)


    # --------------------------------------------------------------------------
    def _run_extract(self):
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)

        print(f">>> Extracting unique values from '{fy.Color.yellow(path_in)}'")
        vgt.VGOperations.extract(path_in)


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
        vgt.VGOperations.rotate(
            path_in, rotate_xy, rotate_yz, rotate_xz
        ).save(path_out)


    # --------------------------------------------------------------------------
    def _run_segment(self):
        path_in    = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out   = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)
        isovalue   = self.main.get_arg_float("isovalue")
        volume_thr = self.main.get_arg_int("volume_thr")

        fmt_in = vg.GridIO.detect_format(path_in)
        if fmt_in != vg.GridFormat.BIN:
            cmd_example = fy.Color.yellow(f"volgrids vgtools convert {path_in} --format BIN")
            self.main.help_and_exit(1,
                f"Segmentation currently only supports 'BIN' format. " +\
                f"Got '{fy.Color.red(fmt_in.name)}' format instead. " +\
                f"Convert the grid first e.g.\n    {cmd_example}"
            )

        fmt_out = vg.GridIO.detect_format(path_out)
        if fmt_out != vg.GridFormat.BIN:
            path_out = Path(str(path_out).replace(path_out.suffix, ".bin"))
            print(
                f"...>>> Warning: output format '{fy.Color.red(fmt_out.name)}' not supported for segmentation.",
                f"Will save as BIN format to '{fy.Color.blue(path_out)}' instead."
            )

        str_iso = fy.Color.green(f"isovalue={isovalue:2.2f}")
        str_thr = fy.Color.cyan(f"volume threshold={volume_thr}")
        print(f">>> Segmenting {fy.Color.yellow(path_in)} with {str_iso} and {str_thr} voxels")
        stdout = vg.Utils.run_bash(
            "_backend/segmentation/run.sh", path_in, path_out, isovalue, volume_thr,
            err_msg = "Grid segmentation failed"
        )
        print(stdout)


    # --------------------------------------------------------------------------
    def _run_average(self):
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        print(f">>> Averaging CMAP file: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.average(path_in).save(path_out)


    # --------------------------------------------------------------------------
    def _run_std_dev(self):
        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        print(f">>> Standard Deviation for CMAP file: {fy.Color.yellow(path_in)}")
        vgt.VGOperations.std_dev(path_in).save(path_out)


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

        print(f">>> Performing '{fy.Color.yellow(command)}' operation on grids: {fy.Color.red(path_in_0)} with {fy.Color.blue(path_in_1)}")
        for grid in vgt.VGOperations.iter_op_binary(
            path_in_0, path_in_1, operation, interpolate_to_common_box
        ):
            grid.save(path_out)


    # --------------------------------------------------------------------------
    def _run_op_unary(self, command: str):
        operation: callable = {
            "abs": vg.Grid.__abs__,
        }[command]

        path_in  = self.main.get_arg_path("path_in",  assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        print(f">>> Performing '{fy.Color.yellow(command)}' operation on grid: {fy.Color.red(path_in)}")

        for grid in vgt.VGOperations.iter_op_unary(path_in, operation):
            grid.save(path_out)


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
        threshold = self.main.get_arg_float("threshold")

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


    # --------------------------------------------------------------------------
    def _get_valid_fmt_arg(self):
        str_fmt = self.main.get_arg_str("format")

        try:
            return vg.GridFormat.from_str(str_fmt)
        except ValueError as e:
            self.main.help_and_exit(1, str(e))


# //////////////////////////////////////////////////////////////////////////////
