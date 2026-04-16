import volgrids as vg
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerSMUtils(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
            "help"   : ("-h", "--help"),
    }


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: volgrids smutils [resids_nonbp|...] [options...]",
            "Available modes:",
            "  resids_nonbp  - Print the set of non-base-paired residue indices in a given RNA structure. A residue is considered non-base-paired if it does not form a canonical base pair (UA, CG) with any other residue. Requires rnapolis",
            "Run 'volgrids smutils [mode] --help' for more details on each mode.",
        )
        if self._has_param_kwds("help") and not self._has_params_pos():
            self._exit_with_help()

        su.OPERATION = self._safe_get_param_pos(0).lower()
        func: callable = self._safe_map_value(su.OPERATION,
            resids_nonbp = self._parse_resids_nonbp,
        )
        func()


    # --------------------------------------------------------------------------
    def _parse_resids_nonbp(self) -> None:
        self._set_help_str(
            "usage: volgrids smutils convert [path/input/grid] [options...]",
            "Available options:",
            "  -h, --help  Show this help message and exit.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        su.PATH_STRUCT = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "Error: 'resids_nonbp' mode requires a path to an input RNA structure file (e.g. rna.pdb)."
            )
        )


# //////////////////////////////////////////////////////////////////////////////
