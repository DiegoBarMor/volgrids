import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifferApp:
    def __init__(self):
        sm.SmifferArgsParser()
        self._apply_custom_config()

        self.ms = sm.SmifferMolecularSystem(sm.PATH_STRUCTURE, sm.PATH_TRAJECTORY)
        self.trimmer: sm.GridTrimmer = None

        str_mode = "PocketSphere" if self.ms.do_ps else "Whole"
        self.timer = vg.Timer(
            f">>> Now processing '{self.ms.molname}' ({vg.USER_MODE}) in '{str_mode}' mode"
        )


    # --------------------------------------------------------------------------
    def run(self):
        if sm.PATH_APBS is None:
            sm.DO_SMIF_APBS = False

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
        trimming_dists = {}
        if sm.DO_SMIF_HYDROPHILIC:
            trimming_dists["small"] = sm.TRIMMING_DIST_SMALL

        if (
            sm.DO_SMIF_STACKING or
            sm.DO_SMIF_HBA or sm.DO_SMIF_HBD or
            sm.DO_SMIF_HYDROPHOBIC or sm.SAVE_TRIMMING_MASK
        ):
            trimming_dists["mid"] = sm.TRIMMING_DIST_MID

        if sm.DO_SMIF_APBS:
            trimming_dists["large"] = sm.TRIMMING_DIST_LARGE

        self.trimmer = sm.GridTrimmer(self.ms, **trimming_dists)

        if sm.SAVE_TRIMMING_MASK:
            mask = self.trimmer.get_mask("mid")
            reverse = vg.Grid.reverse(mask) # save the points that are NOT trimmed
            reverse.save_data(sm.FOLDER_OUT, f"trimming")


        ### Calculate standard SMIF grids
        if sm.DO_SMIF_STACKING:
            self._calc_smif(sm.GridStacking, "mid", "stacking")

        if sm.DO_SMIF_HBA:
            self._calc_smif(sm.GridHBARing, "mid", "hbacceptors")

        if sm.DO_SMIF_HBD:
            self._process_hbdonors("mid")

        if sm.DO_SMIF_HYDROPHOBIC:
            grid_hphob: sm.GridHydrophobic =\
                self._calc_smif(sm.GridHydrophobic, "mid", "hydrophobic")

        if sm.DO_SMIF_HYDROPHILIC:
            grid_hphil: sm.GridHydrophilic =\
                self._calc_smif(sm.GridHydrophilic, "small", "hydrophilic")

        if sm.DO_SMIF_APBS:
            grid_apbs: sm.GridAPBS =\
                self._calc_smif(sm.GridAPBS, "large", "apbs")


        ### Calculate additional grids
        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = grid_hphob - grid_hphil
            grid_hpdiff.save_data(sm.FOLDER_OUT, "hydrodiff")

        if sm.DO_SMIF_LOG_APBS:
            grid_apbs.apply_logabs_transform()
            grid_apbs.save_data(sm.FOLDER_OUT, "apbslog")


    # --------------------------------------------------------------------------
    def _calc_smif(self, cls_grid: type, key_trimming: str, title: str) -> "vg.Grid":
        grid: vg.Grid = cls_grid(self.ms)
        grid.populate_grid()
        self.trimmer.apply_trimming(grid, key_trimming)
        grid.save_data(sm.FOLDER_OUT, title)
        return grid


    # --------------------------------------------------------------------------
    def _process_hbdonors(self, key_trimming: str):
        grid_hbdring = sm.GridHBDRing(self.ms)
        grid_hbdcone = sm.GridHBDCone(self.ms)
        grid_hbdring.populate_grid()
        grid_hbdcone.populate_grid()
        grid_hbd = grid_hbdring + grid_hbdcone
        self.trimmer.apply_trimming(grid_hbd, key_trimming)
        grid_hbd.save_data(sm.FOLDER_OUT, "hbdonors")


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
