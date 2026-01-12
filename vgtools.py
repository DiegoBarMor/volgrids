"""
Script to run Smiffer from the command line.
Installing the volgrids package is not necessary.
"""

import warnings
import volgrids.vgtools as vgt

if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    vgt.AppVGTools.from_cli().run()
