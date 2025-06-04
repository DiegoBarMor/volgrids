import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Smiffer:
    def __init__(self):
        args = vg.args_smiffer()
        self.metadata = vars(args)
        self.metadata["meta"] = args.out / f"{args.pdb.stem}.meta.json"
        self.metadata["do_traj"] = str(args.traj) != '.'
        vg.save_metadata(self.metadata)

        if self.metadata["whole"]:
            self.ms = vg.MSWhole(self.metadata)
            self.Trimmer = vg.TrimmerWhole
        else:
            self.ms = vg.MSPocketSphere(self.metadata)
            self.Trimmer = vg.TrimmerPocketSphere

        self.timer = vg.Timer(
            f">>> Now processing '{args.pdb.stem}' ({'nucleic' if args.rna else 'protein'}) " +\
            f"in '{'Whole' if args.whole else 'PocketSphere'}' mode"
        )


    # --------------------------------------------------------------------------
    def run(self):
        self.timer.start()

        if self.metadata["do_traj"]: # TRAJECTORY MODE
            if not self.metadata["whole"]:
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
            vg.DO_SMIF_STACKING or
            vg.DO_SMIF_HBA or vg.DO_SMIF_HBD or
            vg.DO_SMIF_HYDROPHOBIC or vg.DO_SMIF_APBS
        ):
            trim_large = self.Trimmer(self.ms, vg.TRIMMING_DIST_LARGE)

        if vg.DO_SMIF_HYDROPHILIC:
            trim_small = self.Trimmer(self.ms, vg.TRIMMING_DIST_SMALL)

        ### Calculate standard SMIF grids
        if vg.DO_SMIF_STACKING:
            pg_stacking = vg.GridStacking(self.ms)
            pg_stacking.run(trim_large)

        if vg.DO_SMIF_HBA:
            pg_hba = vg.GridHBAccepts(self.ms)
            pg_hba.run(trim_large)

        if vg.DO_SMIF_HBD:
            pg_hbd = vg.GridHBDonors(self.ms)
            pg_hbd.run(trim_large)

        if vg.DO_SMIF_HYDROPHOBIC:
            pg_hbphob = vg.GridHydrophobic(self.ms)
            pg_hbphob.run(trim_large)

        if vg.DO_SMIF_HYDROPHILIC:
            pg_hphil = vg.GridHydrophilic(self.ms)
            pg_hphil.run(trim_small)

        if vg.DO_SMIF_APBS:
            pg_apbs = vg.GridAPBS(self.ms)
            pg_apbs.run(trim_large)

        ### Calculate additional grids
        if vg.DO_SMIF_LOG_APBS:
            pg_apbs.apply_logabs_transform()

        if vg.DO_SMIF_HYDRODIFF:
            vg.Grid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")


# //////////////////////////////////////////////////////////////////////////////
