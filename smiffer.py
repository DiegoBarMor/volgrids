import volgrids as vg

################################################################################
if __name__ == "__main__":
    ##### Here you can override any of the default parameters.
    ##### Check the volgrids/__init__.py file for other global values (like th SMIFs' gaussian parameters)
    vg.DO_SMIF_STACKING = True
    vg.DO_SMIF_HBA = True
    vg.DO_SMIF_HBD = True
    vg.DO_SMIF_HYDROPHOBIC = True
    vg.DO_SMIF_HYDROPHILIC = True
    vg.DO_SMIF_APBS = True

    vg.DO_SMIF_LOG_APBS = False
    vg.DO_SMIF_HYDRODIFF = False

    vg.DO_TRIMMING_OCCUPANCY = True
    vg.DO_TRIMMING_SPHERE = True # only applies for pocket-sphere mode
    vg.DO_TRIMMING_RNDS = False  # only applies for pocket-sphere mode

    vg.DO_OUTPUT_DX = False
    vg.DO_OUTPUT_MRC = True
    vg.DO_OUTPUT_CMAP = False
    vg.SAVE_CACHED_MASK = False # saves the logical inverse of the trimming mask

    vg.USE_FIXED_DELTAS = True # whether to use fixed dx,dy,dz and le xres,yres,zres change (or the opposite)
    vg.WARNING_GRID_SIZE = 5.0e7 # if the grid would exceed this amount of points, trigger a warning with possibility to abort

    vg.GRID_DX = 0.25 # deltas used for calculations when USE_FIXED_DELTAS = True (resolutions change)
    vg.GRID_DY = 0.25
    vg.GRID_DZ = 0.25

    vg.GRID_XRES = 200 # resolution used for calculations when USE_FIXED_DELTAS = False (deltas change)
    vg.GRID_YRES = 200
    vg.GRID_ZRES = 200


    # --------------------------------------------------------------------------
    meta = vg.SmifferArgsParser()
    if meta.mode in ["prot", "rna"]:
        vg.Smiffer(meta).run()


################################################################################
