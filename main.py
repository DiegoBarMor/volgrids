import volpot as vp

################################################################################
if __name__ == "__main__":
    args = vp.process_args()
    metadata = vars(args)
    metadata["meta"] = args.out / f"{args.pdb.stem}.meta.json"
    vp.save_metadata(metadata)

    print(f">>> Now processing '{args.pdb.stem}' ({'nucleic' if args.rna else 'protein'})...", flush = True)

    Class_ms = vp.MS_Whole if metadata["whole"] else vp.MS_PocketSphere

    ms_trimming_large = Class_ms(metadata, vp.TRIMMING_DIST_LARGE)
    ms_trimming_small = Class_ms(metadata, vp.TRIMMING_DIST_SMALL)
    pg_stacking  = vp.SPG_Stacking   (ms_trimming_large)
    pg_hba       = vp.SPG_HB_Accepts (ms_trimming_large)
    pg_hbd       = vp.SPG_HB_Donors  (ms_trimming_large)
    pg_hbphob    = vp.SPG_Hydrophobic(ms_trimming_large)
    pg_hphil     = vp.SPG_Hydrophilic(ms_trimming_small)
    pg_apbs      = vp.PPG_APBS       (ms_trimming_large)
    pg_hydrodiff = vp.PotentialGrid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")


################################################################################
