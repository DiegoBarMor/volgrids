"""
Script to run Smiffer from the command line.
Installing the volgrids package is not necessary.
"""

import warnings
import volgrids.veins as ve

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    ve.AppVeins.from_cli().run()
