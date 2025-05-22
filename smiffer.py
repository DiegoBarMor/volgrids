import volpot as vp

################################################################################
if __name__ == "__main__":
    ##### Here you can override any of the default parameters.
    ##### Check the volpot/__init__.py file for the list of global values.
    vp.DO_SMIF_APBS = False
    vp.DO_OUTPUT_MRC = False
    vp.DO_OUTPUT_CMAP = True

    vp.Smiffer().run()


################################################################################
