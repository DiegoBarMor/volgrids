import volgrids.vgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifferCalculator:
    def __init__(self, meta: "sm.SmifferArgsParser"):
        self.meta = meta
        self.meta.save_metadata()

        if meta.do_ps:
            self.ms = vg.MSPocketSphere(meta)
            self.Trimmer = sm.TrimmerPocketSphere
            str_mode = "PocketSphere"
        else:
            self.ms = vg.MSWhole(meta)
            self.Trimmer = sm.TrimmerWhole
            str_mode = "Whole"

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
            trim_large = self.Trimmer(self.ms, sm.TRIMMING_DIST_LARGE)

        # if self.meta.do_cavities:
        #     pg_pocket = sm.GridCavities(self.ms)
        #     pg_pocket.radius_pocket = float(self.meta.cavities) # [WIP] improve
        #     pg_pocket.run(trim_large)
        #     return  # Exit after running the pocket finder

        if sm.DO_SMIF_HYDROPHILIC:
            trim_small = self.Trimmer(self.ms, sm.TRIMMING_DIST_SMALL)

        ### Calculate standard SMIF grids
        if sm.DO_SMIF_STACKING:
            pg_stacking = sm.GridStacking(self.ms)
            pg_stacking.run(trim_large)

        if sm.DO_SMIF_HBA:
            pg_hba = sm.GridHBAccepts(self.ms)
            pg_hba.run(trim_large)

        if sm.DO_SMIF_HBD:
            pg_hbd = sm.GridHBDonors(self.ms)
            pg_hbd.run(trim_large)

        if sm.DO_SMIF_HYDROPHOBIC:
            pg_hbphob = sm.GridHydrophobic(self.ms)
            pg_hbphob.run(trim_large)

        if sm.DO_SMIF_HYDROPHILIC:
            pg_hphil = sm.GridHydrophilic(self.ms)
            pg_hphil.run(trim_small)

        if sm.DO_SMIF_APBS:
            pg_apbs = sm.GridAPBS(self.ms)
            pg_apbs.run(trim_large)

        ### Calculate additional grids
        if sm.DO_SMIF_LOG_APBS:
            pg_apbs.apply_logabs_transform()

        if sm.DO_SMIF_HYDRODIFF:
            vg.Grid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")


# //////////////////////////////////////////////////////////////////////////////
