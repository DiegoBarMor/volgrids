"""
Script to run Smiffer from the command line.
Installing the volgrids package is not necessary.
"""

import warnings
import volgrids.smiffer as sm

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    sm.AppSmiffer.from_cli().run()
