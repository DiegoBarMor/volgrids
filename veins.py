"""
Script to run Smiffer from the command line.
Installing the volgrids package is not necessary.
"""

import warnings
import volgrids as vg
import volgrids.veins as ve

if __name__ == "__main__":
    vg.PATH_DEFAULT_CONFIG = vg.resolve_path_resource(__file__, "config_volgrids.ini")
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    ve.AppVeins.from_cli().run()
