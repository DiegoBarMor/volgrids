import volgrids.vgrids as vg
import volgrids.veins as ve

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerVeins(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
        "help"  : ("-h", "--help"),
        "output": ("-o", "--output"),
        "traj"  : ("-t", "--trajectory"),
        "cutoff": ("-c", "--cutoff"),
    }


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: python3 veins.py [mode] [options...]",
            "Available modes:",
            "  energies - Generate grids to visually represent in space the interaction energies of a molecular system.",
            # "  force  - ", # TODO: Implement force mode
            "Run 'python3 veins.py [mode] --help' for more details on each mode.",
            "Running 'python3 veins.py' without a valid mode will display this help message.",
        )
        if self._has_param_kwds("help") and not self._has_params_pos():
            self._exit_with_help(0)

        vg.USER_MODE = self._safe_get_param_pos(0)
        func: callable = self._safe_map_value(vg.USER_MODE,
            energies = self._parse_energies,
            forces   = self._parse_forces,
        )
        func()


    # --------------------------------------------------------------------------
    def _parse_energies(self) -> None:
        self._set_help_str(
            "usage: python3 veins.py energy [path/input/structure.pdb] [path/input/energies.csv] [options...]",
            "Available options:",
            "-h, --help       Show this help message and exit.",
            "-o, --output     Path to the folder where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --trajectory Path to a trajectory file (e.g. XTC) supported by MDAnalysis. In this case, the energies CSV file contains an energy column for each frame. The header of such columns must start with 'frame'.",
            "-c, --cutoff     Energies below this cutoff will be ignored. Default value: 1e-3.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help(0)

        ve.PATH_STRUCTURE = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input structure file provided. Provide a path to the structure file as first positional argument."
            )
        )

        ve.PATH_ENERGIES_CSV = self._safe_path_file_in(
            self._safe_get_param_pos(2,
               err_msg = "No energies CSV file provided. Provide a path to the energies file as second positional argument."
            )
        )

        ve.FOLDER_OUT = self._safe_path_folder_out(
            self._safe_get_param_kwd("output", 0) if self._has_param_kwds("output") \
            else ve.PATH_STRUCTURE.parent
        )

        if self._has_param_kwds("traj"):
            ve.PATH_TRAJECTORY = self._safe_kwd_file_in("traj")

        if self._has_param_kwds("cutoff"):
            try:
                ve.ENERGY_CUTOFF = float(self._safe_get_param_kwd("cutoff", 0))
            except ValueError:
                self._exit_with_help(-1, "The cutoff value must be a valid float number.")


    # --------------------------------------------------------------------------
    def _parse_forces(self) -> None:
        raise NotImplementedError("VEINS Forces mode is not yet implemented.")


# //////////////////////////////////////////////////////////////////////////////
