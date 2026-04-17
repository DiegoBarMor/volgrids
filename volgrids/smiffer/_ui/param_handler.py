from pathlib import Path

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
        "resids": ("-i", "--resids"),
        "pp"    : ("-pp", "--probe-probe"),
    }


    # --------------------------------------------------------------------------
    def assign_globals(self):
        self._set_help_str(
            "usage: volgrids smiffer [prot|rna|ligand] [path.pdb] [options...]",
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
            f"usage: volgrids smiffer {mode} [path.pdb] [options...]",
            "Available options:",
            "-h, --help      Show this help message and exit.",
            "-o, --output    Folder path where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --traj      File path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs      File path to a cached output of APBS for the respective structure file. Prevents the automatic calculation of APBS performed by volgrids' smiffer.",
            "-b, --table     File path to a .chem table file to use for ligand mode, or to override the default macromolecules' tables.",
            "-c, --config    File path to a configuration file with global settings, to override the default settings (e.g. config_volgrids.ini). Alternatively, it can be a list of `configuration=value` keyword pairs for the global settings e.g. (`DO_SMIF_APBS=true DO_SMIF_STACKING=false`).",
            "-s, --sphere    Activate 'pocket sphere' mode by providing the X, Y, Z coordinates (sphere center) and the sphere radius R for a sphere. If not provided, 'whole' mode is assumed.",
            "-i, --resids    File path to a text file containing the residue indices to consider for SMIF calculations (one-based indexing, space separated). Alternatively, a string of space-separated indices can be passed directly as argument. If not provided, all residues will be considered.",
            "-pp, --probe-probe Generate Probe-Probe (PP) fields: spherical accessibility regions (radius 2.0 Å) around interaction sites. Generates hbaPP, hbdPP, stkPP, hpPP fields following overlap-gio2 approach.",
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
        self._handle_params_resids()
        self._handle_params_sphere()
        self._handle_params_pp()
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
        sm.SPHERE = sm.SphereInfo(x_cog, y_cog, z_cog, radius)


    # --------------------------------------------------------------------------
    def _handle_params_configs(self):
        if not self._has_param_kwds("config"): return

        try:
            configs = self._safe_get_param_kwd_list("config")
        except vg.ParamHandler.MissingArgsError:
            ### [WIP] this is not a good solution, it's part of what needs to be refactored when changing the param handler
            available = "\n    " + "\n    ".join(sorted(vg.KNOWN_CONFIGS))
            print(f"Available configuration keys:{available}")
            exit(0)

        for val in configs:
            if '=' in val:
                vg.STR_CUSTOM_CONFIG += f"{val.replace(' ', '\n')}\n"
                continue

            val_as_path = Path(val)
            if not val_as_path.exists():
                self._exit_with_help(self.InvalidPathError, f"The specified config path '{val}' does not exist.")
            if val_as_path.is_dir():
                self._exit_with_help(self.InvalidPathError, f"The specified config path '{val}' is a folder, but a file was expected.")
            vg.PATHS_CUSTOM_CONFIG.append(val_as_path)

    # --------------------------------------------------------------------------
    def _handle_params_resids(self):
        def _handle_path(val: str) -> list[str]:
            possible_path = Path(resids[0])
            if not possible_path.exists():
                return resids
            if possible_path.is_dir():
                self._exit_with_help(self.InvalidPathError, f"The specified resids path '{val}' is a folder, but a file was expected.")
            return possible_path.read_text().strip().split()

        def _assert_resid(resid: str) -> int:
            if not resid.isdigit():
                self._exit_with_help(self.InvalidParamError, f"Invalid residue index '{resid}' provided for --resids option. All indices must be integers.")
            return resid

        if not self._has_param_kwds("resids"): return

        resids = self._safe_get_param_kwd_list("resids")

        if not resids:
            self._exit_with_help(self.InvalidParamError, "No residue indices provided for --resids option.")

        if len(resids) == 1:
            resids = _handle_path(resids[0])

        splitted = (x for resid in resids for x in resid.split())
        sm.CUSTOM_RESIDS = " ".join(_assert_resid(resid) for resid in splitted)

    # --------------------------------------------------------------------------
    def _handle_params_pp(self):
        """Handle probe-probe (PP) field generation flag."""
        if not self._has_param_kwds("pp"):
            return

        # Enable all PP field generation
        sm.DO_SMIF_PP = True
        sm.DO_SMIF_HBA_PP = True
        sm.DO_SMIF_HBD_PP = True
        sm.DO_SMIF_STACKING_PP = True
        sm.DO_SMIF_HYDRO_PP = True

        print(">>> Enabled Probe-Probe (PP) field generation:")
        print("    • hbaPP: HB acceptor spherical accessibility (radius 2.0 Å)")
        print("    • hbdPP: HB donor spherical accessibility (radius 2.0 Å)")
        print("    • stkPP: Stacking spherical accessibility (radius 2.0 Å)")
        print("    • hpPP: Hydrophobic spherical accessibility (radius 2.0 Å)")


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
