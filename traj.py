import volpot as vp

################################################################################
if __name__ == "__main__":
    args = vp.process_args()
    metadata = vars(args)
    metadata["meta"] = args.out / f"{args.pdb.stem}.meta.json"
    vp.save_metadata(metadata)

    print(f">>> Now processing '{args.pdb.stem}' ({'nucleic' if args.rna else 'protein'})...", flush = True)

    if not metadata["whole"]:
        raise NotImplementedError("PocketSphere not implemented yet for trajectory mode.")

    ms = vp.MSWhole(metadata)

    for _ in ms.system.trajectory:
        ms.frame += 1

        print(f">>> Frame {ms.frame}/{len(ms.system.trajectory)}", flush = True)

        trim_large = vp.TrimmerWhole(ms, vp.TRIMMING_DIST_LARGE)
        # trim_small = Trimmer(ms, vp.TRIMMING_DIST_SMALL)

        pg_stacking = vp.GridStacking(ms)
        # pg_hba      = vp.GridHBAccepts(ms)
        # pg_hbd      = vp.GridHBDonors(ms)
        # pg_hbphob   = vp.GridHydrophobic(ms)
        # pg_hphil    = vp.GridHydrophilic(ms)
        # pg_apbs     = vp.GridAPBS(ms)

        pg_stacking.run(trim_large)
        # pg_hba     .run(trim_large)
        # pg_hbd     .run(trim_large)
        # pg_hbphob  .run(trim_large)
        # pg_hphil   .run(trim_small)
        # pg_apbs    .run(trim_large)
        # if vp.DO_LOG_APBS: pg_apbs.apply_logabs_transform()
        # vp.Grid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")


################################################################################
# python3 traj.py -i data/potentials-test/traj/7VKI.pdb -o data/potentials-test/traj -t data/potentials-test/traj/7VKI.xtc -n -w

# coordset #1; vseries play #2
# coordset stop #1; vseries stop #2
