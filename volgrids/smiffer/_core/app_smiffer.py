import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm

try: import freyacli as fy # to display colored text
except ImportError: from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSmiffer(vg.AppSubcommand):
    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain"):
        super().__init__(app_main)
        self.ms: sm.MolSystemSmiffer
        self.trimmer: sm.Trimmer
        self.cavfinder: sm.CavityFinder
        self.timer: vg.Timer
        self.folder_out: Path
        self.path_traj: Path

        mode = self.main.subcommands.pop(0)
        sm.CURRENT_MOLTYPE = sm.MolType.from_str(mode)

        sm.PATH_STRUCT      = self.main.get_arg_path("path_in")
        self.folder_out     = self.main.get_arg_path("folder_out", default = sm.PATH_STRUCT.parent)
        sm.PATH_APBS        = self.main.get_arg_path("path_apbs")
        self.path_traj      = self.main.get_arg_path("path_traj")
        sm.PATH_CHEM_CUSTOM = self.main.get_arg_path("path_chem")

        self.main.assert_file_in(sm.PATH_STRUCT)
        self.main.assert_file_in(sm.PATH_APBS, allow_none = True)
        self.main.assert_file_in(self.path_traj, allow_none = True)
        self.main.assert_file_in(sm.PATH_CHEM_CUSTOM, allow_none = True)
        self.main.assert_dir_out(self.folder_out)

        self._handle_params_configs()
        self._handle_params_resids()
        self._handle_params_sphere()
        self._handle_params_pp()
        self._assert_traj_apbs()
        self._assert_ligand_has_table()


    # --------------------------------------------------------------------------
    def init_smif_parameters(self):
        """Run after calling `AppMain._load_all_configs`."""
        sm.PARAMS_HPHOB = vg.ParamsGaussianUnivariate(
            mu = sm.MU_HYDROPHOBIC, sigma = sm.SIGMA_HYDROPHOBIC,
        )
        sm.PARAMS_HPHIL = vg.ParamsGaussianUnivariate(
            mu = sm.MU_HYDROPHILIC, sigma = sm.SIGMA_HYDROPHILIC,
        )
        sm.PARAMS_HBA = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_HBA, mu_1 = sm.MU_DIST_HBA,
            cov_00 = sm.SIGMA_ANGLE_HBA**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = sm.SIGMA_DIST_HBA**2,
        )
        sm.PARAMS_HBD_FREE = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_HBD_FREE, mu_1 = sm.MU_DIST_HBD_FREE,
            cov_00 = sm.SIGMA_ANGLE_HBD_FREE**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = sm.SIGMA_DIST_HBD_FREE**2,
        )
        sm.PARAMS_HBD_FIXED = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_HBD_FIXED, mu_1 = sm.MU_DIST_HBD_FIXED,
            cov_00 = sm.SIGMA_ANGLE_HBD_FIXED**2, cov_01 = 0,
            cov_10 = 0,  cov_11 = sm.SIGMA_DIST_HBD_FIXED**2,
        )
        sm.PARAMS_STACK = vg.ParamsGaussianBivariate(
            mu_0 = sm.MU_ANGLE_STACKING, mu_1 = sm.MU_DIST_STACKING,
            cov_00 = sm.COV_STACKING_00, cov_01 = sm.COV_STACKING_01,
            cov_10 = sm.COV_STACKING_10, cov_11 = sm.COV_STACKING_11,
        )

        ### square root of the DIST contribution to sm.COV_STACKING,
        sm.SIGMA_DIST_STACKING = np.sqrt(sm.COV_STACKING_11)

        ### these must be initialized after the configs are loaded,
        ### so it's appropriate to place them in this late-init method
        self.ms = sm.MolSystemSmiffer(sm.PATH_STRUCT, self.path_traj)
        self.trimmer = sm.Trimmer.init_infer_dists(self.ms)
        self.cavfinder = sm.CavityFinder()
        self.timer = vg.Timer(
            f">>> SMIFs {sm.CURRENT_MOLTYPE.name:>4} '{fy.Color.yellow(self.ms.molname)}'"+\
            f" in '{'PocketSphere' if self.ms.do_ps else 'Whole'}' mode"
        )


    # --------------------------------------------------------------------------
    def run(self):
        if sm.CURRENT_MOLTYPE.is_ligand():
            sm.DO_SMIF_APBS = False
            print(f"\n...--- ligand: {fy.Color.red('skipping APBS')} SMIF calculation.", end = ' ', flush = True)

        self.timer.start()

        if self.ms.do_traj: # TRAJECTORY MODE
            print()
            for _ in self.ms.system.trajectory:
                self.ms.frame += 1
                timer_frame = vg.Timer(f"...>>> Frame {self.ms.frame}/{len(self.ms.system.trajectory)}")
                timer_frame.start()
                self._process_grids()
                timer_frame.end()
            self._delete_traj_locks()

        else: # SINGLE PDB MODE
            self._process_grids()

        self.timer.end(text = fy.Color.green("SMIFs"), minus = sm.APBS_ELAPSED_TIME)


    # --------------------------------------------------------------------------
    def _process_grids(self):
        ### APBS must be calculated first and split into two parts,
        ### because it can potentially set vg.PQR_CONTENTS_TEMP (used for trimming)
        if sm.DO_SMIF_APBS:
            smif_apbs: sm.SmifAPBS = self._calc_smif(sm.SmifAPBS)
            if vg.PQR_CONTENTS_TEMP:
                new_ms = sm.MolSystemSmiffer.from_pqr_data(vg.PQR_CONTENTS_TEMP)
                sm.MolSystemSmiffer.copy_attributes_except_system(src = self.ms, dst = new_ms)
                self.trimmer.ms = new_ms

        self.trimmer.trim(self.cavfinder)

        if sm.DO_SMIF_APBS:
            smif_apbs.grid.reshape_as_box(self.trimmer.specific_masks["large"].box)
            self._trim_and_save_smif(
                smif_apbs, key_trimming = "large", title = "apbs"
            )
            del self.trimmer.specific_masks["large"]


        ### Calculate standard SMIF grids
        if sm.DO_SMIF_HYDROPHILIC:
            smif_hphil = self._calc_smif(sm.SmifHydrophilic)
            self._trim_and_save_smif(
                smif_hphil, key_trimming = "small", title = "hydrophilic"
            )
            del self.trimmer.specific_masks["small"]

        if sm.DO_SMIF_HYDROPHOBIC:
            smif_hphob = self._calc_smif(sm.SmifHydrophobic)
            self._trim_and_save_smif(
                smif_hphob, key_trimming = "mid", title = "hydrophobic"
            )

        if sm.DO_SMIF_HBA:
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifHBAccepts),
                key_trimming = "mid", title = "hbacceptors"
            )

        if sm.DO_SMIF_HBD:
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifHBDonors),
                key_trimming = "mid", title = "hbdonors"
            )

        if sm.DO_SMIF_STACKING:
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifStacking),
                key_trimming = "mid", title = "stacking"
            )


        ### [TODO] interface calculations will be moved to smutils
        ### Calculate Probe-Probe (PP) fields - spherical accessibility regions
        if sm.DO_SMIF_HBA_PP:
            ### [TODO] probes shouln't be trimmed
            ### at least not occupancy-trimmed (it defeats the purpose of an occupancy probe)
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifHBAcceptsPP),
                key_trimming = "tiny", title = "hbaPP"
            )

        if sm.DO_SMIF_HBD_PP:
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifHBDonorsPP),
                key_trimming = "tiny", title = "hbdPP"
            )

        if sm.DO_SMIF_STACKING_PP:
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifStackingPP),
                key_trimming = "tiny", title = "stkPP"
            )

        if sm.DO_SMIF_HYDRO_PP:
            self._trim_and_save_smif(
                self._calc_smif(sm.SmifHydroPP),
                key_trimming = "tiny", title = "hpPP"
            )


        ### Calculate / store additional grids
        if sm.SAVE_TRIMMING_MASK:
            mask = self.trimmer.get_mask("mid")
            reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
            sm.Smif.save_data_smif(reverse, self.ms, self.folder_out, "trimming")

        if sm.SAVE_CAVITIES and self.cavfinder.has_data():
            sm.Smif.save_data_smif(self.cavfinder.grid, self.ms, self.folder_out, "cavities")

        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = smif_hphob.grid - smif_hphil.grid
            sm.Smif.save_data_smif(grid_hpdiff, self.ms, self.folder_out, "hydrodiff")

        if sm.DO_SMIF_APBS and sm.DO_SMIF_LOG_APBS:
            smif_apbs.apply_logabs_transform()
            sm.Smif.save_data_smif(smif_apbs.grid, self.ms, self.folder_out, "apbslog")

        if not self.ms.do_traj and vg.PQR_CONTENTS_TEMP:
            path_pqr = self.folder_out / f"{self.ms.molname}.pqr"
            path_pqr.write_text(vg.PQR_CONTENTS_TEMP)


    # --------------------------------------------------------------------------
    def _calc_smif(self, cls_smif: type[sm.Smif]) -> "sm.Smif":
        smif: sm.Smif = cls_smif(self.ms)
        smif.populate_grid()
        return smif


    # --------------------------------------------------------------------------
    def _trim_and_save_smif(self, smif: sm.Smif, key_trimming: str, title: str) -> None:
        self.trimmer.mask_grid(smif, key_trimming)
        self.cavfinder.apply_cavities_weighting(smif)
        smif.save_data_smif(smif.grid, self.ms, self.folder_out, title)


    # --------------------------------------------------------------------------
    def _delete_traj_locks(self):
        if self.path_traj.suffix != ".xtc": return

        preffix = str(self.path_traj.parent / f".{self.path_traj.stem}.xtc_offsets")
        Path(f"{preffix}.lock").unlink(missing_ok = True)
        Path(f"{preffix}.npz").unlink(missing_ok = True)


    # --------------------------------------------------------------------------
    def _handle_params_configs(self):
        configs = self.main.get_arg_str("configs")
        if not configs: return

        if isinstance(configs, bool):
            ### this happens when -c is passed without any value to it
            ### `True` is stored instead of any string value
            available = "\n    " + "\n    ".join(sorted(vg.KNOWN_CONFIGS))
            print(f"Available configuration keys:{available}")
            exit(0)

        for val in configs:
            if '=' in val:
                vg.STR_CUSTOM_CONFIG += f"{val.replace(' ', '\n')}\n"
                continue

            val_as_path = Path(val)
            if not val_as_path.exists(): self.main.help_and_exit(1,
                f"The specified config path '{val}' does not exist."
            )
            if val_as_path.is_dir(): self.main.help_and_exit(1,
                f"The specified config path '{val}' is a folder, but a file was expected."
            )
            vg.PATHS_CUSTOM_CONFIG.append(val_as_path)


    # --------------------------------------------------------------------------
    def _handle_params_resids(self):
        def _handle_possible_path(val: str) -> str:
            possible_path = Path(resids[0])
            if not possible_path.exists(): return resids
            if possible_path.is_dir(): self.main.help_and_exit(1,
                f"The specified resids path '{val}' is a folder, but a file was expected."
            )
            return possible_path.read_text().strip()

        def _assert_resid(resid: str) -> int:
            if not resid.isdigit(): self.main.help_and_exit(1,
                f"Invalid residue index '{resid}' provided for --resids option. All indices must be integers."
            )
            return resid

        resids = self.main.get_arg_str("resids")
        if not resids: return

        if not resids: self.main.help_and_exit(1,
            "No residue indices provided for --resids option."
        )

        if len(resids) == 1: # [NOTE] when passing a path, only one input file is currrently supported
            resids = _handle_possible_path(resids[0])

        splitted = (x for resid in resids for x in resid.split())
        sm.CUSTOM_RESIDS = " ".join(_assert_resid(resid) for resid in splitted)


    # --------------------------------------------------------------------------
    def _handle_params_sphere(self):
        sphere = self.main.get_arg_list_float("sphere")
        if not sphere: return

        sm.SPHERE = sm.SphereInfo(*sphere)


    # --------------------------------------------------------------------------
    def _handle_params_pp(self):
        """Handle probe-probe (PP) field generation flag."""
        if not self.main.get_arg_bool("probe_probe"): return

        ### [TODO] controlling the PP field generation both as a config and through a flag could be confusing to use
        # Enable all PP field generation
        sm.DO_SMIF_HBA_PP = True
        sm.DO_SMIF_HBD_PP = True
        sm.DO_SMIF_STACKING_PP = True
        sm.DO_SMIF_HYDRO_PP = True

        ### [TODO] it's fine for proof-of-concept but having 5 extra rows printed every time could clutter the user's output
        ### could be replaced with a little 'probes' or 'interfaces' mode instead e.g. `SMIFs PROT 'name' in 'probes' mode`
        ### that approach would need some changes on how that print is handled though
        print(">>> Enabled Probe-Probe (PP) field generation:")
        print("    • hbaPP: HB acceptor spherical accessibility (radius 2.0 Å)")
        print("    • hbdPP: HB donor spherical accessibility (radius 2.0 Å)")
        print("    • stkPP: Stacking spherical accessibility (radius 2.0 Å)")
        print("    • hpPP: Hydrophobic spherical accessibility (radius 2.0 Å)")


    # --------------------------------------------------------------------------
    def _assert_traj_apbs(self):
        if not sm.DO_SMIF_APBS: return
        if (self.path_traj is None) or (sm.PATH_APBS is None): return
        self.main.help_and_exit(1,
            f"The APBS output '{sm.PATH_APBS}' was provided. However, "+
            "trajectory mode is enabled, so this file would be ambiguous. "+
            "Please either disable trajectory mode or remove the APBS file input. "
        )


    # --------------------------------------------------------------------------
    def _assert_ligand_has_table(self):
        if not sm.CURRENT_MOLTYPE.is_ligand(): return
        if sm.PATH_CHEM_CUSTOM is not None: return
        self.main.help_and_exit(1,
            "No table file provided for ligand mode. Use -b or --table to specify the path to the .chem table file."
        )


# //////////////////////////////////////////////////////////////////////////////
