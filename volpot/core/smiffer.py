import volpot as vp

# //////////////////////////////////////////////////////////////////////////////
class Smiffer:
    def __init__(self):
        args = vp.args_smiffer()
        self.metadata = vars(args)
        self.metadata["meta"] = args.out / f"{args.pdb.stem}.meta.json"
        self.metadata["do_traj"] = str(args.traj) != '.'
        vp.save_metadata(self.metadata)

        if self.metadata["whole"]:
            self.ms = vp.MSWhole(self.metadata)
            self.Trimmer = vp.TrimmerWhole
        else:
            self.ms = vp.MSPocketSphere(self.metadata)
            self.Trimmer = vp.TrimmerPocketSphere

        self.timer = vp.Timer(
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
                timer_frame = vp.Timer(f"...>>> Frame {self.ms.frame}/{len(self.ms.system.trajectory)}")
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
            vp.DO_SMIF_STACKING or
            vp.DO_SMIF_HBA or vp.DO_SMIF_HBD or
            vp.DO_SMIF_HYDROPHOBIC or vp.DO_SMIF_APBS
        ):
            trim_large = self.Trimmer(self.ms, vp.TRIMMING_DIST_LARGE)

        if vp.DO_SMIF_HYDROPHILIC:
            trim_small = self.Trimmer(self.ms, vp.TRIMMING_DIST_SMALL)

        ### Calculate standard SMIF grids
        if vp.DO_SMIF_STACKING:
            pg_stacking = vp.GridStacking(self.ms)
            pg_stacking.run(trim_large)

        if vp.DO_SMIF_HBA:
            pg_hba = vp.GridHBAccepts(self.ms)
            pg_hba.run(trim_large)

        if vp.DO_SMIF_HBD:
            pg_hbd = vp.GridHBDonors(self.ms)
            pg_hbd.run(trim_large)

        if vp.DO_SMIF_HYDROPHOBIC:
            pg_hbphob = vp.GridHydrophobic(self.ms)
            pg_hbphob.run(trim_large)

        if vp.DO_SMIF_HYDROPHILIC:
            pg_hphil = vp.GridHydrophilic(self.ms)
            pg_hphil.run(trim_small)

        if vp.DO_SMIF_APBS:
            pg_apbs = vp.GridAPBS(self.ms)
            pg_apbs.run(trim_large)

        ### Calculate additional grids
        if vp.DO_SMIF_LOG_APBS:
            pg_apbs.apply_logabs_transform()

        if vp.DO_SMIF_HYDRODIFF:
            vp.Grid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")


# //////////////////////////////////////////////////////////////////////////////
