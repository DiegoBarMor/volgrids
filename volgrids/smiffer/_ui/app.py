import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class AppSmiffer(vg.App):
    CONFIG_MODULES = (vg, sm)
    _CLASS_PARAM_HANDLER = sm.ParamHandlerSmiffer

    _CLASS_TRIMMER = sm.Trimmer
    _CLASS_MOL_SYSTEM = sm.MolSystemSmiffer

    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._init_globals()

        self.ms: sm.MolSystemSmiffer = self._CLASS_MOL_SYSTEM(sm.PATH_STRUCT, sm.PATH_TRAJ)
        self.trimmer: sm.Trimmer = self._CLASS_TRIMMER.init_infer_dists(self.ms)
        self.cavfinder: sm.CavityFinder = sm.CavityFinder()
        self.timer = vg.Timer(
            f">>> SMIFs {sm.CURRENT_MOLTYPE.name:>4} '{self.ms.molname}'"+\
            f" in '{'PocketSphere' if self.ms.do_ps else 'Whole'}' mode"
        )


    # --------------------------------------------------------------------------
    def run(self):
        if sm.CURRENT_MOLTYPE.is_ligand():
            sm.DO_SMIF_APBS = False
            print("\n...--- ligand: skipping APBS SMIF calculation.", end = ' ', flush = True)

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

        self.timer.end(text = "SMIFs", minus = sm.APBS_ELAPSED_TIME)


    # --------------------------------------------------------------------------
    def _init_globals(self):
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


        ### Calculate Probe-Probe (PP) fields - spherical accessibility regions
        if sm.DO_SMIF_HBA_PP:
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
            sm.Smif.save_data_smif(reverse, self.ms, sm.FOLDER_OUT, "trimming")

        if sm.SAVE_CAVITIES and self.cavfinder.has_data():
            sm.Smif.save_data_smif(self.cavfinder.grid, self.ms, sm.FOLDER_OUT, "cavities")

        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = smif_hphob.grid - smif_hphil.grid
            sm.Smif.save_data_smif(grid_hpdiff, self.ms, sm.FOLDER_OUT, "hydrodiff")

        if sm.DO_SMIF_APBS and sm.DO_SMIF_LOG_APBS:
            smif_apbs.apply_logabs_transform()
            sm.Smif.save_data_smif(smif_apbs.grid, self.ms, sm.FOLDER_OUT, "apbslog")

        if not self.ms.do_traj and vg.PQR_CONTENTS_TEMP:
            path_pqr = sm.FOLDER_OUT / f"{self.ms.molname}.pqr"
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
        smif.save_data_smif(smif.grid, self.ms, sm.FOLDER_OUT, title)


    # --------------------------------------------------------------------------
    def _delete_traj_locks(self):
        if sm.PATH_TRAJ.suffix != ".xtc": return

        preffix = str(sm.PATH_TRAJ.parent / f".{sm.PATH_TRAJ.stem}.xtc_offsets")
        Path(f"{preffix}.lock").unlink(missing_ok = True)
        Path(f"{preffix}.npz").unlink(missing_ok = True)


# //////////////////////////////////////////////////////////////////////////////
