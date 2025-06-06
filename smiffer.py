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

    vg.WARNING_GRID_SIZE = 5.0e7 # if the grid would exceed this amount of points, trigger a warning with possibility to abort

    vg.Smiffer().run()


################################################################################
