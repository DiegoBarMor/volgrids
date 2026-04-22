import volgrids as vg
import volgrids.smutils as su

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


    # --------------------------------------------------------------------------
    def _run_occupancy(self) -> None:
        app = su.AppOccupancy(self.main)
        app.init_smif_parameters()
        app.run()


    # --------------------------------------------------------------------------
    def _run_resids_nonbp(self) -> None:
        path_in = self.main.get_arg_path("path_in")
        self.main.assert_file_in(path_in)
        su.SMOperations.print_resids_nonbp(path_in)


# //////////////////////////////////////////////////////////////////////////////
