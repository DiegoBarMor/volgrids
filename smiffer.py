import volgrids as vg

################################################################################
if __name__ == "__main__":
    ##### Here you can override any of the default parameters.
    ##### Check the volgrids/__init__.py file for the list of global values.
    vg.DO_SMIF_APBS = False
    vg.DO_OUTPUT_MRC = False
    vg.DO_OUTPUT_CMAP = True

    vg.Smiffer().run()


################################################################################
