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

        self.timer.end()


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
        self.trimmer.trim()

        if sm.SAVE_TRIMMING_MASK:
            mask = self.trimmer.get_mask("mid")
            reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
            reverse.save_data(sm.FOLDER_OUT, f"trimming")

        ### Calculate standard SMIF grids
        if sm.DO_SMIF_STACKING:
            self._calc_smif(sm.SmifStacking, "mid", "stacking")

        if sm.DO_SMIF_HBA:
            self._calc_smif(sm.SmifHBAccepts, "mid", "hbacceptors")

        if sm.DO_SMIF_HBD:
            self._calc_smif(sm.SmifHBDonors, "mid", "hbdonors")

        if sm.DO_SMIF_HYDROPHOBIC:
            grid_hphob: sm.SmifHydrophobic =\
                self._calc_smif(sm.SmifHydrophobic, "mid", "hydrophobic")

        if sm.DO_SMIF_HYDROPHILIC:
            grid_hphil: sm.SmifHydrophilic =\
                self._calc_smif(sm.SmifHydrophilic, "small", "hydrophilic")

        if sm.DO_SMIF_APBS:
            grid_apbs: sm.SmifAPBS =\
                self._calc_smif(sm.SmifAPBS, "large", "apbs")


        ### Calculate additional grids
        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = grid_hphob - grid_hphil
            grid_hpdiff.save_data(sm.FOLDER_OUT, "hydrodiff")

        if sm.DO_SMIF_LOG_APBS:
            grid_apbs.apply_logabs_transform()
            grid_apbs.save_data(sm.FOLDER_OUT, "apbslog")


    # --------------------------------------------------------------------------
    def _calc_smif(self, cls_grid: type[sm.Smif], key_trimming: str, title: str) -> "vg.Grid":
        grid: sm.Smif = cls_grid(self.ms)
        grid.populate_grid()
        self.trimmer.mask_grid(grid, key_trimming)
        grid.save_data(sm.FOLDER_OUT, title)
        return grid


# //////////////////////////////////////////////////////////////////////////////
