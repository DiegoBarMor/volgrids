import volpot as vp

################################################################################
if __name__ == "__main__":
    args = vp.process_args()
    metadata = vars(args)
    metadata["meta"] = args.out / f"{args.pdb.stem}.meta.json"
    vp.save_metadata(metadata)

    timer = vp.Timer(
        f">>> Now processing '{args.pdb.stem}' ({'nucleic' if args.rna else 'protein'}) " +\
        f"in '{'Whole' if metadata['whole'] else 'PocketSphere'}' mode"
    )

    if metadata["whole"]:
        ms = vp.MSWhole(metadata)
        Trimmer = vp.TrimmerWhole
    else:
        ms = vp.MSPocketSphere(metadata)
        Trimmer = vp.TrimmerPocketSphere

    trim_large = Trimmer(ms, vp.TRIMMING_DIST_LARGE)
    trim_small = Trimmer(ms, vp.TRIMMING_DIST_SMALL)

    pg_stacking = vp.GridStacking(ms)
    pg_hba      = vp.GridHBAccepts(ms)
    pg_hbd      = vp.GridHBDonors(ms)
    pg_hbphob   = vp.GridHydrophobic(ms)
    pg_hphil    = vp.GridHydrophilic(ms)
    pg_apbs     = vp.GridAPBS(ms)

    pg_stacking.run(trim_large)
    pg_hba     .run(trim_large)
    pg_hbd     .run(trim_large)
    pg_hbphob  .run(trim_large)
    pg_hphil   .run(trim_small)
    pg_apbs    .run(trim_large)
    if vp.DO_LOG_APBS: pg_apbs.apply_logabs_transform()
    vp.Grid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")

    timer.end()


################################################################################
