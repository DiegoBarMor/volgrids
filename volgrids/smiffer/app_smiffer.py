import numpy as np
from pathlib import Path
import multiprocessing as mp
from collections import deque

import volgrids as vg
import volgrids.smiffer as sm
from volgrids._vendors import freyacli as fy

### Set by the parent before spawning the trajectory MP pool; inherited by workers via fork.
_WORKER_APP: "AppSmiffer" = None


# ------------------------------------------------------------------------------
def _worker_init():
    """Re-create the MolSystem in each worker so file handles aren't shared across forks."""
    _WORKER_APP.ms = sm.MolSystemSmiffer(sm.PATH_STRUCT, _WORKER_APP.path_traj)


# ------------------------------------------------------------------------------
def _worker_process_frame(frame_idx: int) -> tuple[int, float]:
    return _WORKER_APP._process_frame(frame_idx)


# //////////////////////////////////////////////////////////////////////////////
class AppSmiffer(vg.AppSubcommand):
    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain", str_mode: str = "SMIFs"):
        super().__init__(app_main)
        self.str_mode = str_mode # just for printing

        self.ms: sm.MolSystemSmiffer
        self.trimmer: sm.Trimmer
        self.cavfinder: sm.CavityFinder
        self.timer: vg.Timer
        self.folder_out: Path
        self.path_traj: Path
        self.nproc: int

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
        self.ms = sm.MolSystemSmiffer(sm.PATH_STRUCT, self.path_traj)
        self.trimmer = sm.Trimmer.init_infer_dists(self.ms)
        self.cavfinder = sm.CavityFinder()
        self.timer = vg.Timer(
            f">>> {sm.CURRENT_MOLTYPE.name:<4} {'Sphere' if self.ms.do_ps else 'Whole'} "+\
            f"{fy.Color.magenta(self.str_mode)} for '{fy.Color.yellow(self.ms.molname)}'"
        )


    # --------------------------------------------------------------------------
    def run(self):
        def _end(is_traj: bool):
            self.timer.end(text = fy.Color.green("volgrids"), minus = sm.APBS_ELAPSED_TIME)
            if is_traj: self._delete_traj_locks()

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
            self._run_traj_parallel(n_frames)
            return _end(is_traj = True)

        ### 1.b) TRAJECTORY MODE (single process)
        for i in range(n_frames):
            self._process_frame(i)
        return _end(is_traj = True)


    # --------------------------------------------------------------------------
    def _process_frame(self, frame_idx: int) -> None:
        """Set per-frame state and run `_process_grids`."""
        self.ms.system.trajectory[frame_idx]
        self.ms.frame = frame_idx + 1
        self.trimmer = sm.Trimmer.init_infer_dists(self.ms)
        self.cavfinder = sm.CavityFinder()

        n_frames = len(self.ms.system.trajectory)
        timer = vg.Timer(f"...>>> Frame {self.ms.frame}/{n_frames}")
        timer.start()
        self._process_grids()
        timer.end()


    # --------------------------------------------------------------------------
    def _run_traj_parallel(self, n_frames: int) -> None:
        global _WORKER_APP

        ### Pre-clear stale CMAP outputs once in the parent. Workers must not race on this.
        if vg.REMOVE_OLD_CMAP_OUTPUT:
            for path in self.folder_out.glob(f"{self.ms.molname}.*.cmap"):
                path.unlink()
            vg.REMOVE_OLD_CMAP_OUTPUT = False

        ### "fork" so children inherit module state (configs, _WORKER_APP, etc.) without pickling
        ctx = mp.get_context("fork")
        vg.MP_CMAP_LOCK = ctx.Lock()
        _WORKER_APP = self

        nproc = min(self.nproc, n_frames)
        try:
            with ctx.Pool(nproc, initializer = _worker_init) as pool:
                deque(pool.imap_unordered(_worker_process_frame, range(n_frames)), maxlen = 0)
        finally:
            vg.MP_CMAP_LOCK = None
            _WORKER_APP = None


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
