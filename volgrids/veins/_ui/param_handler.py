import volgrids as vg
import volgrids.veins as ve

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerVeins(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
        "help"  : ("-h", "--help"),
        "output": ("-o", "--output"),
        "traj"  : ("-t", "--trajectory"),
    }


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: volgrids veins [mode] [options...]",
            "Available modes:",
            "  energy - Generate grids for the spatial interaction energies of a molecular system.",
            "  force   - Generate grids for the spatial interaction forces of a molecular system.",
            "Run 'volgrids veins [mode] --help' for more details on each mode.",
        )
        if self._has_param_kwds("help") and not self._has_params_pos():
            self._exit_with_help()

        ve.MODE = self._safe_get_param_pos(0).lower()
        func: callable = self._safe_map_value(ve.MODE,
            energy = self._parse_energy,
            force  = self._parse_force,
        )
        func()


    # --------------------------------------------------------------------------
    def _parse_energy(self) -> None:
        raise NotImplementedError("[WIP]")

        self._set_help_str(
            "usage: volgrids veins energy [path/input/structure.pdb] [path/input/energies.csv] [options...]",
            "Available options:",
            "-h, --help       Show this help message and exit.",
            "-o, --output     Path to the folder where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --trajectory Enable trajectory mode. Use this flag to provide a frame number that will be used for the grid's key in the output CMAP file.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        ve.PATH_CSV_IN = self._safe_path_file_in(
            self._safe_get_param_pos(2,
               err_msg = "No energies CSV file provided. Provide a path to the energies file as second positional argument."
            )
        )

        ve.FOLDER_OUT = self._safe_kwd_folder_out("output", default = ve.PATH_CSV_IN.parent)
        ve.TRAJ_FRAME = self._safe_kwd_int("traj", default = None)


    # --------------------------------------------------------------------------
    def _parse_force(self) -> None:
        self._set_help_str(
            "usage: volgrids veins force [path/input/structure.pdb] [path/input/forces.csv] [options...]",
            "Available options:",
            "-h, --help       Show this help message and exit.",
            "-o, --output     Path to the folder where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --trajectory Enable trajectory mode. Use this flag to provide a frame number that will be used for the grid's key in the output CMAP file.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        ve.PATH_CSV_IN = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No forces CSV file provided. Provide a path to the forces file as second positional argument."
            )
        )

        ve.FOLDER_OUT = self._safe_kwd_folder_out("output", default = ve.PATH_CSV_IN.parent)
        ve.TRAJ_FRAME = self._safe_kwd_int("traj", default = None)
        if ve.TRAJ_FRAME is not None:
            raise NotImplementedError("[WIP]")


# //////////////////////////////////////////////////////////////////////////////
