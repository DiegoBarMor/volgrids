import volgrids as vg
import volgrids.vgtools as vgt

DEFAULT_COMPARISON_THRESHOLD = 1e-5 # [TODO] this should be a config

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.AppSubcommand):
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


    # --------------------------------------------------------------------------
    def _run_convert(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["out_dx", "out_mrc", "out_ccp4", "out_cmap"],
            allow_none = True,
        )
        path_in       = self.main.get_arg_path("path_in")
        path_out_dx   = self.main.get_arg_path("out_dx")
        path_out_mrc  = self.main.get_arg_path("out_mrc")
        path_out_ccp4 = self.main.get_arg_path("out_ccp4")
        path_out_cmap = self.main.get_arg_path("out_cmap")

        def _convert(path_out, fmt_out: vg.GridFormat):
            if path_out is None: return
            print(f">>> Converting {path_in} file to {fmt_out.name}: {path_out}")
            vgt.VGOperations.convert(path_in, path_out, fmt_out)

        _convert(path_out_dx,   vg.GridFormat.DX)
        _convert(path_out_mrc,  vg.GridFormat.MRC)
        _convert(path_out_ccp4, vg.GridFormat.CCP4)
        _convert(path_out_cmap, vg.GridFormat.CMAP)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        self.main.assert_paths(
            keys_file_in = ["paths_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        paths_in = self.main.get_arg_list_path("paths_in")
        path_out = self.main.get_arg_path("path_out")

        print(f">>> Packing {len(paths_in)} grids into '{path_out}'")
        vgt.VGOperations.pack(paths_in, path_out)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        self.main.assert_paths(keys_file_in = ["path_in"], allow_none = False)
        path_in = self.main.get_arg_path("path_in")

        if self.main.get_arg_value("folder_out") is None:
            self.main.set_arg_value("folder_out", path_in.parent)

        self.main.assert_paths(keys_dir_out = ["folder_out"], allow_none = True)
        path_out = self.main.get_arg_path("folder_out")

        print(f">>> Unpacking '{path_in}' into '{path_out}'")
        vgt.VGOperations.unpack(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_fix_cmap(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")

        print(f">>> Fixing CMAP file: {path_in}")
        vgt.VGOperations.fix_cmap(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_average(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            keys_file_out = ["path_out"],
            allow_none = False,
        )
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")

        print(f">>> Averaging CMAP file: {path_in}")
        vgt.VGOperations.average(path_in, path_out)


    # --------------------------------------------------------------------------
    def _run_summary(self):
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            allow_none = False,
        )
        path_in = self.main.get_arg_path("path_in")

        print(f">>> Grid summary: {path_in}")
        vgt.VGOperations.summary(path_in)


    # --------------------------------------------------------------------------
    def _run_compare(self):
        self.main.assert_paths(
            keys_file_in = ["path_0", "path_1"],
            allow_none = False,
        )
        path_in_0 = self.main.get_arg_path("path_0")
        path_in_1 = self.main.get_arg_path("path_1")

        threshold = self.main.get_arg_float("threshold")
        if threshold is None:
            threshold = DEFAULT_COMPARISON_THRESHOLD

        print(f">>> Comparing grids: {path_in_0} vs {path_in_1} (threshold={threshold:2.2e})")
        result = vgt.VGOperations.compare(path_in_0, path_in_1, threshold)

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
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")

        rotate_yz = self.main.get_arg_float("x")
        rotate_xz = self.main.get_arg_float("y")
        rotate_xy = self.main.get_arg_float("z")

        print(f">>> Rotating grid: {path_in} by {rotate_xy}° (xy), {rotate_yz}° (yz), {rotate_xz}° (xz)")
        vgt.VGOperations.rotate(path_in, path_out, rotate_xy, rotate_yz, rotate_xz)


# //////////////////////////////////////////////////////////////////////////////
