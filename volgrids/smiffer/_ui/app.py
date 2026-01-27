import numpy as np

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
        self.timer = vg.Timer(
            f">>> Now processing {sm.CURRENT_MOLTYPE.name:>4} '{self.ms.molname}'"+\
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

        else: # SINGLE PDB MODE
            self._process_grids()

        self.timer.end(text = "SMIFS", minus = sm.APBS_ELAPSED_TIME)


    # --------------------------------------------------------------------------
    def _import_config_dependencies(self):
        return {"np": np, "vg": vg}


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
                self.trimmer.ms = sm.MolSystemSmiffer.from_pqr_data(vg.PQR_CONTENTS_TEMP)


        self.trimmer.trim()

        if sm.DO_SMIF_APBS:
            smif_apbs.reshape_as_other(self.trimmer.specific_masks["large"])
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


        ### Calculate additional grids
        if sm.SAVE_TRIMMING_MASK:
            mask = self.trimmer.get_mask("mid")
            reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
            reverse.save_data(sm.FOLDER_OUT, f"trimming")

        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = smif_hphob - smif_hphil
            grid_hpdiff.save_data(sm.FOLDER_OUT, "hydrodiff")

        if sm.DO_SMIF_LOG_APBS:
            smif_apbs.apply_logabs_transform()
            smif_apbs.save_data(sm.FOLDER_OUT, "apbslog")


    # --------------------------------------------------------------------------
    def _calc_smif(self, cls_smif: type[sm.Smif]) -> "sm.Smif":
        smif: sm.Smif = cls_smif(self.ms)
        smif.populate_grid()
        return smif


    # --------------------------------------------------------------------------
    def _trim_and_save_smif(self, smif: sm.Smif, key_trimming: str, title: str) -> None:
        self.trimmer.mask_grid(smif, key_trimming)
        smif.save_data(sm.FOLDER_OUT, title)


# //////////////////////////////////////////////////////////////////////////////
