import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSmiffer(vg.AppSubcommand):
    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain", str_mode: str = "SMIFs"):
        super().__init__(app_main)
        self.str_mode = str_mode # just for printing

        ##### initialized once in __init__
        self.ms: sm.MolSystem
        self.timer: vg.Timer
        self.folder_out: Path
        self.path_traj: Path
        self.nproc: int

        #### set in every call of _process_grids
        self.trimmer: sm.Trimmer
        self.cavfinder: sm.CavityFinder
        self.grid_smif: vg.Grid

        mode = self.main.subcommands.pop(0)
        sm.CURRENT_MOLTYPE = sm.MolType.from_str(mode)

        sm.PATH_STRUCT      = self.main.get_arg_path("path_in",   assertion = fy.PathAssertion.FILE_IN)
        sm.PATH_APBS        = self.main.get_arg_path("path_apbs", assertion = fy.PathAssertion.FILE_IN)
        self.path_traj      = self.main.get_arg_path("path_traj", assertion = fy.PathAssertion.FILE_IN)
        sm.PATH_CHEM_CUSTOM = self.main.get_arg_path("path_chem", assertion = fy.PathAssertion.FILE_IN)
        self.folder_out     = self.main.get_arg_path("folder_out",
            default = sm.PATH_STRUCT.parent, assertion = fy.PathAssertion.DIR_OUT
        )
        self.nproc          = max(1, self.main.get_arg_int("nproc", default = 1))

        self._handle_params_configs()
        self._handle_params_resids()
        self._handle_params_sphere()
        self._assert_traj_apbs()
        self._assert_ligand_has_table()

        app_main.load_configs(vg, sm)

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
        ### so it's appropriate to place them at the end of init
        self.ms = sm.MolSystem(sm.PATH_STRUCT, self.path_traj)
        self.timer = vg.Timer(
            f">>> {sm.CURRENT_MOLTYPE.name:<4} {'Sphere' if self.ms.do_ps else 'Whole'} "+\
            f"{fy.Color.magenta(self.str_mode)} for '{fy.Color.yellow(self.ms.molname)}'"
        )


    # --------------------------------------------------------------------------
    def run(self):
        def _end(is_traj: bool):
            self.timer.end(text = fy.Color.green("volgrids"), minus = sm.APBS_ELAPSED_TIME)
            if is_traj:
                self._delete_traj_locks()
            elif vg.TMP_APBS_CONTENT_PQR:
                path_pqr = self.folder_out / f"{self.ms.molname}.pqr"
                path_pqr.write_text(vg.TMP_APBS_CONTENT_PQR)


        if sm.CURRENT_MOLTYPE.is_ligand():
            sm.DO_SMIF_APBS = False
            print(f"\n...--- ligand: {fy.Color.red('skipping APBS')} SMIF calculation.", end = ' ', flush = True)

        self.timer.start()

        ##### 0) SINGLE PDB MODE
        if not self.ms.do_traj:
            self._process_grids()
            return _end(is_traj = False)
        #####

        print()
        n_frames = len(self.ms.system.trajectory)

        ### 1.a) TRAJECTORY MODE (multiprocessing)
        if self.nproc > 1:
            sm.TrajMultiprocess(self).run(n_frames)
            return _end(is_traj = True)

        ### 1.b) TRAJECTORY MODE (single process)
        for i in range(n_frames):
            self.process_frame(i)
        return _end(is_traj = True)


    # --------------------------------------------------------------------------
    def process_frame(self, frame_idx: int) -> None:
        """Set per-frame state and run `_process_grids`."""
        self.ms.system.trajectory[frame_idx]
        self.ms.frame = frame_idx + 1

        n_frames = len(self.ms.system.trajectory)
        timer = vg.Timer(f"...>>> Frame {self.ms.frame}/{n_frames}")
        timer.start()
        self._process_grids()
        timer.end()


    # --------------------------------------------------------------------------
    def _process_grids(self):
        ms = self._current_mol_system()
        self.trimmer = sm.Trimmer.init_infer_dists(ms)
        self.cavfinder = sm.CavityFinder()

        self.trimmer.trim(self.cavfinder)

        if sm.DO_SMIF_APBS:
            ### [TODO] refactor into taking directly box (instead of initializing an empty grid_apbs)
            ### this is implemented like this to ave a consistent signature for populate_grid among different smifs
            ### it's also the only reason why populate_grid returns a grid
            grid_apbs = vg.Grid(ms.box, init_grid = False)
            grid_apbs = sm.SmifAPBS(ms).populate_grid(grid_apbs)
            grid_apbs.reshape_as_box(self.trimmer.specific_masks["large"].box) # [TODO] needed?
            self._trim_and_weight_grid(grid_apbs, key_trimming = "large")
            sm.Smif.save_data(grid_apbs, ms, self.folder_out, "apbs")
            del grid_apbs
            del self.trimmer.specific_masks["large"]

        self.grid_smif = vg.Grid(ms.box)

        ### Calculate standard SMIF grids
        if sm.DO_SMIF_HYDROPHILIC:
            sm.SmifHydrophilic(ms).populate_grid(self.grid_smif)
            self._trim_and_weight_grid(self.grid_smif, key_trimming = "small")
            sm.Smif.save_data(self.grid_smif, ms, self.folder_out, "hydrophilic")
            del self.trimmer.specific_masks["small"]

        if sm.DO_SMIF_HYDROPHOBIC:
            sm.SmifHydrophobic(ms).populate_grid(self.grid_smif)
            self._trim_and_weight_grid(self.grid_smif, key_trimming = "mid")
            sm.Smif.save_data(self.grid_smif, ms, self.folder_out, "hydrophobic")

        if sm.DO_SMIF_HBA:
            sm.SmifHBAccepts(ms).populate_grid(self.grid_smif)
            self._trim_and_weight_grid(self.grid_smif, key_trimming = "mid")
            sm.Smif.save_data(self.grid_smif, ms, self.folder_out, "hbacceptors")

        if sm.DO_SMIF_HBD:
            sm.SmifHBDonors(ms).populate_grid(self.grid_smif)
            self._trim_and_weight_grid(self.grid_smif, key_trimming = "mid")
            sm.Smif.save_data(self.grid_smif, ms, self.folder_out, "hbdonors")

        if sm.DO_SMIF_STACKING:
            sm.SmifStacking(ms).populate_grid(self.grid_smif)
            self._trim_and_weight_grid(self.grid_smif, key_trimming = "mid")
            sm.Smif.save_data(self.grid_smif, ms, self.folder_out, "stacking")


        ### Calculate / store additional grids
        if sm.SAVE_TRIMMING_MASK:
            mask = self.trimmer.get_mask("mid")
            reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
            sm.Smif.save_data(reverse, ms, self.folder_out, "trimming")

        if sm.SAVE_CAVITIES and self.cavfinder.has_data():
            sm.Smif.save_data(self.cavfinder.grid, ms, self.folder_out, "cavities")

        del self.grid_smif


    # --------------------------------------------------------------------------
    def _current_mol_system(self) -> "sm.MolSystem":
        if not sm.DO_SMIF_APBS: return self.ms

        ### When dealing with APBS, a first step of PQR generation should always be performed.
        ### PQR can add hydrogen atoms, which should be reflected in the occupancy trimming
        ### operation or be used by HBond SMIFs.
        sm.SmifAPBS(self.ms).gen_pqr()
        obj = sm.MolSystem.from_pqr_data(vg.TMP_APBS_CONTENT_PQR)
        sm.MolSystem.copy_attributes_except_system(src = self.ms, dst = obj)
        return obj


    # --------------------------------------------------------------------------
    def _trim_and_weight_grid(self, grid: vg.Grid, key_trimming: str) -> None:
        self.trimmer.mask_grid(grid, key_trimming)
        self.cavfinder.apply_cavities_weighting(grid)


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
                vg.STR_CUSTOM_CONFIG += val.replace(' ', '\n') + "\n"
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
        def _assert_residue(residue: str) -> int:
            splitted = residue.split('.')
            if len(splitted) != 2: self.main.help_and_exit(1,
                f"Invalid residue '{residue}' provided for --residues option." +\
                "The residues must be in the format \"chain_id.resid\" (e.g. \"A.3 A.4 A.5 B.10\")"
            )
            return residue

        residues = self.main.get_arg_str("residues")
        if not residues: return

        splitted = (x for res in residues for x in res.split())
        sm.CUSTOM_RESIDUES = " ".join(_assert_residue(res) for res in splitted)

        if not sm.CUSTOM_RESIDUES: print(
            fy.Color.red("WARNING: ") + "No valid residues provided for --residues option."
        )


    # --------------------------------------------------------------------------
    def _handle_params_sphere(self):
        sphere = self.main.get_arg_float("sphere", is_list = True)
        if not sphere: return

        sm.SPHERE = sm.SphereInfo(*sphere)


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
