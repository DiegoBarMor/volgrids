import volgrids.vgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class ParamHandlerSmiffer(vg.ParamHandler):
    _EXPECTED_CLI_FLAGS = {
        "help"  : ("-h", "--help"),
        "output": ("-o", "--output"),
        "traj"  : ("-t", "--traj"),
        "apbs"  : ("-a", "--apbs"),
        "pocket": ("-rxyz", "-ps", "--pocket-sphere"),
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
            self._exit_with_help(0)

        vg.USER_MODE = self._safe_get_param_pos(0)
        sm.CURRENT_MOLTYPE = self._safe_map_value(vg.USER_MODE,
            prot = sm.MolType.PROT,
            rna = sm.MolType.RNA,
            ligand = sm.MolType.LIGAND,
        )


        self._set_help_str(
            f"usage: python3 smiffer.py {vg.USER_MODE} [path/input/struct.pdb] [options...]",
            "Available options:",
            "-h, --help                       Show this help message and exit.",
            "-o, --output                     Folder path where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --traj                       File path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs                       File path to the output of APBS for the respective structure file (this must be done before). An OpenDX file is expected.",
            "-b, --table                      File path to a .chem table file to use for ligand mode, or to override the default macromolecules' tables.",
            "-c, --config                     File path to a configuration file with global settings, to override the default settings from config.ini.",
            "-rxyz", "-ps", "--pocket-sphere  Activate 'pocket sphere' mode by providing the sphere radius and the X, Y, Z coordinates for its center. If not provided, 'whole' mode is assumed.",
        )
        if self._has_param_kwds("help"):
            self._exit_with_help(0)

        sm.PATH_STRUCTURE = self._safe_path_file_in(
            self._safe_get_param_pos(1,
               err_msg = "No input structure file provided. Provide a path to the structure file as first positional argument."
            )
        )

        sm.FOLDER_OUT = self._safe_path_folder_out(
            self._safe_get_param_kwd("output", 0) if self._has_param_kwds("output") \
            else sm.PATH_STRUCTURE.parent
        )

        if self._has_param_kwds("traj"):
            sm.PATH_TRAJECTORY = self._safe_kwd_file_in("traj")

        if self._has_param_kwds("apbs"):
            sm.PATH_APBS = self._safe_kwd_file_in("apbs")

        if self._has_param_kwds("table"):
            sm.PATH_TABLE = self._safe_kwd_file_in("table")
        elif sm.CURRENT_MOLTYPE == sm.MolType.LIGAND:
            self._exit_with_help(-1, "No table file provided for ligand mode. Use -b or --table to specify the path to the .chem table file.")

        if self._has_param_kwds("config"):
            sm.PATH_CONFIG = self._safe_kwd_file_in("config")

        if self._has_param_kwds("pocket"):
            params_pocket = self._params_kwd["pocket"]
            try:
                radius = float(self._safe_idx(params_pocket, 0, "Missing pocket sphere radius."))
                x_cog  = float(self._safe_idx(params_pocket, 1, "Missing pocket sphere center X coordinate."))
                y_cog  = float(self._safe_idx(params_pocket, 2, "Missing pocket sphere center Y coordinate."))
                z_cog  = float(self._safe_idx(params_pocket, 3, "Missing pocket sphere center Z coordinate."))
            except ValueError:
                self._exit_with_help(-1, "Pocket sphere options must be numeric values.")
            sm.PS_INFO = (radius, x_cog, y_cog, z_cog)


# //////////////////////////////////////////////////////////////////////////////
