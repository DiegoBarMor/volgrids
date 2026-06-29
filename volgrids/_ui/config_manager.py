import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class ConfigManager:
    def __init__(self):
        ##################### VOLGRIDS
        self.out_format: str = "MRC"
        self.out_cmap_compression: int = 9
        self.out_warning_npoints: float = 1.0e99
        self.out_overwrite_ok: bool = True

        self.box_padding: int = 5
        self.box_tight_traj: bool = False
        self.box_force_equilateral: bool = False
        self.box_from_deltas: bool = True

        self.box_dx: float = 0.25
        self.box_dy: float = 0.25
        self.box_dz: float = 0.25

        self.box_xres: int = 200
        self.box_yres: int = 200
        self.box_zres: int = 200

        ##################### SMIFFER
        self.smif_stk:   bool = True
        self.smif_hba:   bool = True
        self.smif_hbd:   bool = True
        self.smif_apbs:  bool = True
        self.smif_hphob: bool = True
        self.smif_hphil: bool = True
        self.smif_use_hydrogens: bool = True
        self.smif_hb_only_nbase: bool = False

        self.trim_sphere:    bool = True
        self.trim_occupancy: bool = True
        self.trim_rnds:      bool = True
        self.trim_faraway:   bool = True
        self.trim_cavities:  bool = False
        self.trim_save: bool = False

        self.cav_threshold: int = 3
        self.cav_npasses: int = 2
        self.cav_weight: float = 0.0
        self.cav_save: bool = False

        self.trim_occ_dist_short: float = 2.5
        self.trim_occ_dist_mid:   float = 3.0
        self.trim_occ_dist_long:  float = 3.5

        self.trim_rnds_max_dist:    float = float("inf")
        self.trim_rnds_cube_radius: int = 4

        self.trim_faraway_dist: float = 7.0

        self.param_stk_scale: float = 3.5
        self.param_hb_scale:  float = 3.5

        self.param_hphob_dist_mu:    float = 4.4
        self.param_hbhob_dist_sigma: float = 0.3

        self.param_hphil_dist_mu:    float = 3.0
        self.param_hphil_dist_sigma: float = 0.15

        self.param_hba_angle_mu:    float = 129.9
        self.param_hba_dist_mu:     float = 2.93
        self.param_hba_angle_sigma: float = 20.0
        self.param_hba_dist_sigma:  float = 0.21

        self.param_hbd_free_angle_mu:    float = 109.0
        self.param_hbd_free_dist_mu:     float = 2.93
        self.param_hbd_free_angle_sigma: float = 20.0
        self.param_hbd_free_dist_sigma:  float = 0.21

        self.param_hbd_fixed_angle_mu:    float = 180.0
        self.param_hbd_fixed_dist_mu:     float = 2.93
        self.param_hbd_fixed_angle_sigma: float = 30.0
        self.param_hbd_fixed_dist_sigma:  float = 0.21

        self.param_stk_angle_mu: float = 29.97767535
        self.param_stk_dist_mu:  float = 4.1876158
        self.param_stk_cov00:    float = 169.9862228
        self.param_stk_cov01:    float = 6.62318852
        self.param_stk_cov10:    float = 6.62318852
        self.param_stk_cov11:    float = 0.37123882

        self.misc_kernel_gaussian_sigmas: int = 4
        self.misc_logapbs_min: int = -2
        self.misc_logapbs_max: int = 3

        ##################### SMUTILS
        self.og_stk_radius: float = 2.0
        self.og_hba_radius: float = 2.0
        self.og_hbd_radius: float = 2.0
        # self.og_apbs_radius: float = 2.0
        self.og_hphob_radius: float = 2.0
        # self.og_hphil_radius: float = 2.0

        self.debug_chemtable_ligand: bool = False


    # --------------------------------------------------------------------------
    def update_configs_from_ini(self, ini: vg.ParserIni) -> None:
        headers = ini.headers()
        if len(headers) > 1:
            raise ValueError(
                f"Configuration file must have no headers, but the following were found: {headers[1:]}. "+\
                "Remove all headers and keep only key-value pairs (or comments) in the file."
            )
        assert headers[0] == '', f"Unexpected INI header with value '{headers[0]}'. Check ParserIni."

        for key, value in ini.iter_splitted_lines(''):
            vg.CFG.register_config(key, value)


    # --------------------------------------------------------------------------
    def register_config(self, key: str, str_value: str) -> None:
        def parse_str_value():
            ### INTEGERS
            if str_value.isdigit():
                return int(str_value)

            ### FLOATS
            try: return float(str_value)
            except ValueError: pass

            ### BOOLEANS
            if str_value.lower() in ["true", "false"]:
                return str_value.lower() == "true"

            ### STRINGS
            return str_value.strip('"').strip("'")

        if key.upper() not in vg.KNOWN_CONFIGS:
            raise ValueError(f"Unknown configuration: {key}.")

        self.__dict__[key.lower()] = parse_str_value()


    # --------------------------------------------------------------------------
    def display_help(self):
        available = "\n    " + "\n    ".join(sorted(vg.KNOWN_CONFIGS))
        print(f"Available configuration keys:{available}")
        exit(0)


# //////////////////////////////////////////////////////////////////////////////
