import volgrids as vg
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class AppSMUtils(vg.AppSubcommand):
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.func_operation: callable = None


    # --------------------------------------------------------------------------
    def assign_globals(self):
        su.OPERATION = self.main.subcommands.pop(0)
        self.func_operation = {
            "resids_nonbp" : self._run_resids_nonbp,
        }[su.OPERATION]


    # --------------------------------------------------------------------------
    def run(self):
        self.func_operation()


    # --------------------------------------------------------------------------
    def _run_resids_nonbp(self) -> None:
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            allow_none = False,
        )
        su.PATH_STRUCT = self.main.get_arg_path("path_in")

        su.SMOperations.print_resids_nonbp(su.PATH_STRUCT)


# //////////////////////////////////////////////////////////////////////////////
