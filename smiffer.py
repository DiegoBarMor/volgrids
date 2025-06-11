import volgrids.vgrids as vg
import volgrids.smiffer as sm

################################################################################
if __name__ == "__main__":
    ##### Here you can override any of the default parameters.
    ##### Check the volgrids/__init__.py file for other global values (like th SMIFs' gaussian parameters)
    sm.DO_SMIF_STACKING = True
    sm.DO_SMIF_HBA = True
    sm.DO_SMIF_HBD = True
    sm.DO_SMIF_HYDROPHOBIC = True
    sm.DO_SMIF_HYDROPHILIC = True
    sm.DO_SMIF_APBS = True

    sm.DO_SMIF_LOG_APBS = False
    sm.DO_SMIF_HYDRODIFF = False

    sm.DO_TRIMMING_OCCUPANCY = True
    sm.DO_TRIMMING_SPHERE = True # only applies for pocket-sphere mode
    sm.DO_TRIMMING_RNDS = False  # only applies for pocket-sphere mode

    vg.DO_OUTPUT_DX = False # Human-readable format (heavy). Tested with VMD, Chimera, ChimeraX, UnityMol.
    vg.DO_OUTPUT_MRC = False # Binary format (light). Tested with VMD, Chimera, ChimeraX.
    vg.DO_OUTPUT_CMAP = True # Compressed binary format (very light). Tested with ChimeraX.
    sm.SAVE_CACHED_MASK = False # saves the logical inverse of the trimming mask

    vg.USE_FIXED_DELTAS = True # whether to use fixed dx,dy,dz and le xres,yres,zres change (or the opposite)
    vg.WARNING_GRID_SIZE = 5.0e7 # if the grid would exceed this amount of points, trigger a warning with possibility to abort

    vg.GRID_DX = 0.25 # deltas used for calculations when USE_FIXED_DELTAS = True (resolutions change)
    vg.GRID_DY = 0.25
    vg.GRID_DZ = 0.25

    vg.GRID_XRES = 200 # resolution used for calculations when USE_FIXED_DELTAS = False (deltas change)
    vg.GRID_YRES = 200
    vg.GRID_ZRES = 200


    # --------------------------------------------------------------------------
    meta = sm.SmifferArgsParser()

    for name, value in meta.debug_vars.items():
        if hasattr(vg, name):
            setattr(vg, name, value)
        elif hasattr(sm, name):
            setattr(sm, name, value)
        else:
            print(f"Warning: {name} is not a valid volgrids variable. Skipping.")
            continue

    if meta.mode in ["prot", "rna"]:
        sm.SmifferCalculator(meta).run()
    else:
        vg.VGTools(meta).run()


################################################################################
