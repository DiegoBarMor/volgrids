import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerSmiffer(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
        "help"  : ("-h", "--help"),
        "output": ("-o", "--output"),
        "traj"  : ("-t", "--traj"),
        "apbs"  : ("-a", "--apbs"),
        "sphere": ("-s", "--sphere"),
        "table" : ("-b", "--table"),
        "config": ("-c", "--config"),
    }


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: volgrids smiffer [prot|rna|convert|pack|unpack] [options...]",
            "Available modes:",
            "  prot     - Calculate SMIFs for protein structures.",
            "  rna      - Calculate SMIFs for RNA structures.",
            "  ligand   - Calculate SMIFs for ligand structures. A .chem table must be provided.",
            "Run 'volgrids smiffer [mode] --help' for more details on each mode.",
        )
        if self._has_param_kwds("help") and not self._has_params_pos():
            self._exit_with_help()

        mode = self._safe_get_param_pos(0)
        sm.CURRENT_MOLTYPE = self._safe_map_value(mode.lower(),
            prot = sm.MolType.PROT,
            rna = sm.MolType.RNA,
            ligand = sm.MolType.LIGAND,
        )

        self._set_help_str(
            f"usage: volgrids smiffer {mode} [path/input/struct.pdb] [options...]",
            "Available options:",
            "-h, --help        Show this help message and exit.",
            "-o, --output      Folder path where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --traj        File path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs        File path to a cached output of APBS for the respective structure file. Prevents the automatic calculation of APBS performed by volgrids' smiffer.",
            "-b, --table       File path to a .chem table file to use for ligand mode, or to override the default macromolecules' tables.",
            "-c, --config      File path to a configuration file with global settings, to override the default settings (e.g. config_volgrids.ini).",
            "-s, --sphere      Activate 'pocket sphere' mode by providing the X, Y, Z coordinates (sphere center) and the sphere radius R for a sphere. If not provided, 'whole' mode is assumed.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        sm.PATH_STRUCT = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input structure file provided. Provide a path to the structure file as first positional argument."
            )
        )

        sm.FOLDER_OUT = self._safe_kwd_folder_out("output", default = sm.PATH_STRUCT.parent)
        sm.PATH_APBS  = self._safe_kwd_file_in("apbs")
        sm.PATH_TRAJ  = self._safe_kwd_file_in("traj")
        sm.PATH_TABLE = self._safe_kwd_file_in("table")

        self._handle_params_configs()
        self._handle_params_sphere()
        self._assert_traj_apbs()
        self._assert_ligand_has_table()


    # --------------------------------------------------------------------------
    def _handle_params_sphere(self):
        if not self._has_param_kwds("sphere"): return
        params_sphere = self._params_kwd["sphere"]
        try:
            x_cog  = float(self._safe_idx(params_sphere, 0, "Missing sphere center X coordinate."))
            y_cog  = float(self._safe_idx(params_sphere, 1, "Missing sphere center Y coordinate."))
            z_cog  = float(self._safe_idx(params_sphere, 2, "Missing sphere center Z coordinate."))
            radius = float(self._safe_idx(params_sphere, 3, "Missing sphere radius."))
        except ValueError:
            self._exit_with_help(self.InvalidParamError, "Sphere options must be numeric values.")
        sm.SPHERE_INFO = (x_cog, y_cog, z_cog, radius)


    # --------------------------------------------------------------------------
    def _handle_params_configs(self):
        path_config = self._safe_kwd_file_in("config") # [WIP]
        if path_config is None: return
        vg.PATHS_CUSTOM_CONFIG = [path_config]


    # --------------------------------------------------------------------------
    def _assert_traj_apbs(self):
        if not sm.DO_SMIF_APBS: return
        if (sm.PATH_TRAJ is None) or (sm.PATH_APBS is None): return
        raise self.InvalidParamError(
            f"The APBS output '{sm.PATH_APBS}' was provided. However, "+
            "trajectory mode is enabled, so this file would be ambiguous. "+
            "Please either disable trajectory mode or remove the APBS file input. "
        )


    # --------------------------------------------------------------------------
    def _assert_ligand_has_table(self):
        if not sm.CURRENT_MOLTYPE.is_ligand(): return
        if self._has_param_kwds("table"): return
        self._exit_with_help(
            self.MissingParamError,
            "No table file provided for ligand mode. Use -b or --table to specify the path to the .chem table file."
        )


# //////////////////////////////////////////////////////////////////////////////
