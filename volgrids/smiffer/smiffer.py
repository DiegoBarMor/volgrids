import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifferApp:
    def __init__(self):
        sm.SmifferArgsParser()
        self._apply_custom_config()

        self.ms = sm.SmifferMolecularSystem(sm.PATH_STRUCTURE, sm.PATH_TRAJECTORY)

        str_mode = "PocketSphere" if self.ms.do_ps else "Whole"
        self.timer = vg.Timer(
            f">>> Now processing '{self.ms.molname}' ({vg.USER_MODE}) in '{str_mode}' mode"
        )


    # --------------------------------------------------------------------------
    def run(self):
        self.timer.start()

        if self.ms.do_traj: # TRAJECTORY MODE
            if self.ms.do_ps:
                raise NotImplementedError("PocketSphere not implemented yet for trajectory mode. Use -w flag")

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
    def _process_grids(self):
        ### Only trim if needed
        if (
            sm.DO_SMIF_STACKING or
            sm.DO_SMIF_HBA or sm.DO_SMIF_HBD or
            sm.DO_SMIF_HYDROPHOBIC or sm.DO_SMIF_APBS or
            sm.SAVE_CACHED_MASK
        ):
            trim_large = sm.GridTrimmer(self.ms, sm.TRIMMING_DIST_LARGE)

        if sm.DO_SMIF_HYDROPHILIC:
            trim_small = sm.GridTrimmer(self.ms, sm.TRIMMING_DIST_SMALL)

        ### Calculate standard SMIF grids
        if sm.DO_SMIF_STACKING:
            self._calc_smif(sm.GridStacking(self.ms), trim_large, "stacking")

        if sm.DO_SMIF_HBA:
            self._calc_smif(sm.GridHBARing(self.ms), trim_large, "hbacceptors")

        if sm.DO_SMIF_HBD:
            grid_hbdring = sm.GridHBDRing(self.ms)
            grid_hbdcone = sm.GridHBDCone(self.ms)
            grid_hbdring.populate_grid()
            grid_hbdcone.populate_grid()
            grid_hbd = grid_hbdring + grid_hbdcone
            # grid_hbd = grid_hbdcone
            trim_large.apply_trimming(grid_hbd)
            grid_hbd.save_data(sm.FOLDER_OUT, "hbdonors")

        if sm.DO_SMIF_HYDROPHOBIC:
            grid_hphob = sm.GridHydrophobic(self.ms)
            self._calc_smif(grid_hphob, trim_large, "hydrophobic")

        if sm.DO_SMIF_HYDROPHILIC:
            grid_hphil = sm.GridHydrophilic(self.ms)
            self._calc_smif(grid_hphil, trim_small, "hydrophilic")

        if sm.DO_SMIF_APBS:
            self._process_apbs(trim_large)

        ### Calculate additional grids
        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = vg.Grid.substract(grid_hphob, grid_hphil)
            grid_hpdiff.save_data(sm.FOLDER_OUT, "hydrodiff")


    # --------------------------------------------------------------------------
    def _calc_smif(self, grid: "vg.Grid", trimmer: "sm.GridTrimmer", title: str):
        grid.populate_grid()
        trimmer.apply_trimming(grid)
        grid.save_data(sm.FOLDER_OUT, title)


    # --------------------------------------------------------------------------
    def _process_apbs(self, trimmer: "sm.GridTrimmer"):
        if sm.PATH_APBS is None: return

        grid_apbs = sm.GridAPBS(self.ms)
        self._calc_smif(grid_apbs, trimmer, "apbs")

        if not sm.DO_SMIF_LOG_APBS: return

        grid_apbs.apply_logabs_transform()
        grid_apbs.save_data(sm.FOLDER_OUT, "apbslog")


    # --------------------------------------------------------------------------
    def _apply_custom_config(self) -> None:
        if sm.PATH_CONFIG is None: return

        config = vg.ConfigParser(sm.PATH_CONFIG)
        if config.has("VOLGRIDS"): config.apply_config(
            key = "VOLGRIDS", scope = vg.__dict__,
            valid_configs = set(vg.__annotations__.keys()),
            all_configs_mandatory = False
        )
        if config.has("SMIFFER"): config.apply_config(
            key = "SMIFFER", scope = sm.__dict__,
            valid_configs = set(sm.__annotations__.keys()),
            all_configs_mandatory = False
        )


# //////////////////////////////////////////////////////////////////////////////
