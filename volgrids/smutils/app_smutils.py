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
        if operation == "resids_nonbp": return self._run_resids_nonbp()
        if operation == "plot_dist"   : return self._run_plot_dist()


    # --------------------------------------------------------------------------
    def _run_occupancy(self) -> None:
        su.AppOccupancy(self.main).run()


    # --------------------------------------------------------------------------
    def _run_resids_nonbp(self) -> None:
        path_in = self.main.get_arg_path("path_in")
        self.main.assert_file_in(path_in)
        resids_nonbp = su.RNAResids.get_resids_nonbp(path_in)
        print(resids_nonbp)


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


# //////////////////////////////////////////////////////////////////////////////
