import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSMUtils(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.func_operation: callable = None


    # --------------------------------------------------------------------------
    def run(self):
        operation = self.main.subcommands.pop(0)
        if operation == "occupancy": return self._run_occupancy()
        if operation == "pwoverlap": return self._run_pwoverlap()
        if operation == "res_nobp" : return self._run_res_nobp()
        if operation == "res_nostk": return self._run_res_nostk()
        if operation == "histogram": return self._run_histogram()
        if operation == "log_apbs" : return self._run_log_apbs()
        raise ValueError(f"Unknown operation: {operation}")


    # --------------------------------------------------------------------------
    def _run_occupancy(self) -> None:
        su.AppOccupancy(self.main).run()


    # --------------------------------------------------------------------------
    def _run_pwoverlap(self) -> None:
        su.AppPwOverlap(self.main).run()


    # --------------------------------------------------------------------------
    def _run_res_nobp(self) -> None:
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        print(su.ResiduesNucleic.get_residues_nobp(path_in))


    # --------------------------------------------------------------------------
    def _run_res_nostk(self) -> None:
        path_in = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        print(su.ResiduesNucleic.get_residues_nostk(path_in))


    # --------------------------------------------------------------------------
    def _run_histogram(self) -> None:
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT, allow_none = True)
        key      = self.main.get_arg_str("key")

        print(f">>> Voxel distribution: {fy.Color.yellow(path_in)}")
        su.Histogram.plot(path_in, path_out, key)


    # --------------------------------------------------------------------------
    def _run_log_apbs(self) -> None:
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_out = self.main.get_arg_path("path_out", assertion = fy.PathAssertion.FILE_OUT)

        grid = vg.GridIO.read_auto(path_in)
        sm.SmifAPBS.apply_logabs_transform(grid)
        vg.GridIO.write_auto(grid, path_out)


# //////////////////////////////////////////////////////////////////////////////
