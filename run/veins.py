"""
Script to run Smiffer from the command line.
Installing the volgrids package is not necessary.
Please run from the root of the repository, e.g.:
    python run/smiffer.py --help
"""

import warnings
from _mods.allow_root_imports import *
import volgrids.veins as ve

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    ve.AppVeins.from_cli().run()
