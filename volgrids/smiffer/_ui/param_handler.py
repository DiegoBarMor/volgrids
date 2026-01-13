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
            "usage: python3 smiffer.py [prot|rna|convert|pack|unpack] [options...]",
            "Available modes:",
            "  prot     - Calculate SMIFs for protein structures.",
            "  rna      - Calculate SMIFs for RNA structures.",
            "  ligand   - Calculate SMIFs for ligand structures. A .chem table must be provided.",
            "Run 'python3 smiffer.py [mode] --help' for more details on each mode.",
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
            f"usage: python3 smiffer.py {mode} [path/input/struct.pdb] [options...]",
            "Available options:",
            "-h, --help        Show this help message and exit.",
            "-o, --output      Folder path where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --traj        File path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs        File path to the output of APBS for the respective structure file (this must be done before). An OpenDX file is expected.",
            "-b, --table       File path to a .chem table file to use for ligand mode, or to override the default macromolecules' tables.",
            "-c, --config      File path to a configuration file with global settings, to override the default settings (e.g. config_volgrids.ini).",
            "-s, --sphere      Activate 'pocket sphere' mode by providing the X, Y, Z coordinates (sphere center) and the sphere radius R for a sphere. If not provided, 'whole' mode is assumed.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help()

        if sm.CURRENT_MOLTYPE.is_ligand() and not self._has_param_kwds("table"):
            self._exit_with_help("No table file provided for ligand mode. Use -b or --table to specify the path to the .chem table file.")

        sm.PATH_STRUCTURE = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input structure file provided. Provide a path to the structure file as first positional argument."
            )
        )

        sm.FOLDER_OUT = self._safe_kwd_folder_out("output", default = sm.PATH_STRUCTURE.parent)
        sm.PATH_TRAJECTORY    = self._safe_kwd_file_in("traj")
        sm.PATH_APBS          = self._safe_kwd_file_in("apbs")
        sm.PATH_TABLE         = self._safe_kwd_file_in("table")
        vg.PATH_CUSTOM_CONFIG = self._safe_kwd_file_in("config")

        if self._has_param_kwds("sphere"):
            params_sphere = self._params_kwd["sphere"]
            try:
                x_cog  = float(self._safe_idx(params_sphere, 0, "Missing sphere center X coordinate."))
                y_cog  = float(self._safe_idx(params_sphere, 1, "Missing sphere center Y coordinate."))
                z_cog  = float(self._safe_idx(params_sphere, 2, "Missing sphere center Z coordinate."))
                radius = float(self._safe_idx(params_sphere, 3, "Missing sphere radius."))
            except ValueError:
                self._exit_with_help("Sphere options must be numeric values.")
            sm.SPHERE_INFO = (x_cog, y_cog, z_cog, radius)


# //////////////////////////////////////////////////////////////////////////////
