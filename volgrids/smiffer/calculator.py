import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifferCalculator:
    def __init__(self, meta: "sm.SmifferArgsParser"):
        self.meta = meta
        self.meta.save_metadata()
        self.ms = sm.SmifferMolecularSystem(meta)

        str_mode = "PocketSphere" if meta.do_ps else "Whole"
        self.timer = vg.Timer(
            f">>> Now processing '{meta.name}' ({meta.mode}) in '{str_mode}' mode"
        )


    # --------------------------------------------------------------------------
    def run(self):
        self.timer.start()

        if self.meta.do_traj: # TRAJECTORY MODE
            if self.meta.do_ps:
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
        # do_cavities = self.meta.do_cavities

        ### Only trim if needed
        if (
            sm.DO_SMIF_STACKING or
            sm.DO_SMIF_HBA or sm.DO_SMIF_HBD or
            sm.DO_SMIF_HYDROPHOBIC or sm.DO_SMIF_APBS or
            sm.SAVE_CACHED_MASK # or do_cavities
        ):
            trim_large = sm.GridTrimmer(self.ms, sm.TRIMMING_DIST_LARGE)

        # if self.meta.do_cavities:
        #     pg_pocket = sm.GridCavities(self.ms)
        #     pg_pocket.radius_pocket = float(self.meta.cavities) # [WIP] improve
        #     pg_pocket.run(trim_large)
        #     return  # Exit after running the pocket finder

        if sm.DO_SMIF_HYDROPHILIC:
            trim_small = sm.GridTrimmer(self.ms, sm.TRIMMING_DIST_SMALL)

        ### Calculate standard SMIF grids
        if sm.DO_SMIF_STACKING:
            self._calc_smif(sm.GridStacking(self.ms), trim_large)

        if sm.DO_SMIF_HBA:
            self._calc_smif(sm.GridHBAccepts(self.ms), trim_large)

        if sm.DO_SMIF_HBD:
            self._calc_smif(sm.GridHBDonors(self.ms), trim_large)

        if sm.DO_SMIF_HYDROPHOBIC:
            grid_hphob = sm.GridHydrophobic(self.ms)
            self._calc_smif(grid_hphob, trim_large)

        if sm.DO_SMIF_HYDROPHILIC:
            grid_hphil = sm.GridHydrophilic(self.ms)
            self._calc_smif(grid_hphil, trim_small)

        if sm.DO_SMIF_APBS:
            self._process_apbs(trim_large)

        ### Calculate additional grids
        if sm.DO_SMIF_HYDROPHOBIC and sm.DO_SMIF_HYDROPHILIC and sm.DO_SMIF_HYDRODIFF:
            grid_hpdiff = vg.Grid.substract(grid_hphob, grid_hphil)
            grid_hpdiff.save_data(override_prefix = "hydrodiff")


    # --------------------------------------------------------------------------
    def _calc_smif(self, grid: "vg.Grid", trimmer: "sm.GridTrimmer"):
        grid.populate_grid()
        trimmer.apply_trimming(grid)
        grid.save_data()


    # --------------------------------------------------------------------------
    def _process_apbs(self, trimmer: "sm.GridTrimmer"):
        if self.meta.path_apbs is None: return

        grid_apbs = sm.GridAPBS(self.ms)
        self._calc_smif(grid_apbs, trimmer)

        if not sm.DO_SMIF_LOG_APBS: return

        grid_apbs.apply_logabs_transform()
        grid_apbs.save_data(override_prefix = "apbslog")


# //////////////////////////////////////////////////////////////////////////////
