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
        if operation == "resids_nonbp": return self._run_resids_nonbp()


    # --------------------------------------------------------------------------
    def _run_resids_nonbp(self) -> None:
        self.main.assert_paths(
            keys_file_in = ["path_in"],
            allow_none = False,
        )
        su.SMOperations.print_resids_nonbp(
            self.main.get_arg_path("path_in")
        )


# //////////////////////////////////////////////////////////////////////////////
