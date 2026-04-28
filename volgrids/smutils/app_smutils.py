import volgrids as vg
import volgrids.smutils as su

try: import freyacli as fy
except ImportError: from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSMUtils(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.func_operation: callable = None


    # --------------------------------------------------------------------------
    def run(self):
        operation = self.main.subcommands.pop(0)
        if operation == "occupancy"   : return self._run_occupancy()
        if operation == "resids_nobp" : return self._run_resids_nobp()
        if operation == "resids_nostk": return self._run_resids_nostk()
        if operation == "histogram"   : return self._run_plot_dist()
        raise ValueError(f"Unknown operation: {operation}")


    # --------------------------------------------------------------------------
    def _run_occupancy(self) -> None:
        su.AppOccupancy(self.main).run()


    # --------------------------------------------------------------------------
    def _run_resids_nobp(self) -> None:
        path_in = self.main.get_arg_path("path_in")
        self.main.assert_file_in(path_in)
        resids = su.RNAResids.get_resids_nobp(path_in)
        print(resids)


    # --------------------------------------------------------------------------
    def _run_resids_nostk(self) -> None:
        path_in = self.main.get_arg_path("path_in")
        self.main.assert_file_in(path_in)
        resids = su.RNAResids.get_resids_nostk(path_in)
        print(resids)


    # --------------------------------------------------------------------------
    def _run_plot_dist(self) -> None:
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")
        key      = self.main.get_arg_str("key")

        self.main.assert_file_in(path_in)
        if path_out is not None:
            self.main.assert_file_out(path_out)

        print(f">>> Voxel distribution: {fy.Color.yellow(path_in)}")
        su.Histogram.plot(path_in, path_out, key)


    # --------------------------------------------------------------------------
    def _run_plot_3d(self) -> None:
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")
        key      = self.main.get_arg_str("key")
        stride   = self.main.get_arg_int("stride", default = 2)

        self.main.assert_file_in(path_in)
        if path_out is not None:
            self.main.assert_file_out(path_out)

        if isinstance(key, list): key = key[0] if key else None

        print(f">>> 3D volume: {fy.Color.yellow(path_in)}")
        su.DistPlot.plot_3d(path_in, path_out, key, stride)


    # --------------------------------------------------------------------------
    def _run_plot_echarts(self) -> None:
        path_in  = self.main.get_arg_path("path_in")
        path_out = self.main.get_arg_path("path_out")
        key      = self.main.get_arg_str("key")
        stride   = self.main.get_arg_int("stride", default = 1)

        self.main.assert_file_in(path_in)
        if path_out is not None:
            self.main.assert_file_out(path_out)

        if isinstance(key, list): key = key[0] if key else None

        print(f">>> ECharts 3D scatter: {fy.Color.yellow(path_in)}")
        su.DistPlot.plot_echarts(path_in, path_out, key, stride)


# //////////////////////////////////////////////////////////////////////////////
