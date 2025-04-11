import argparse
from argparse import RawDescriptionHelpFormatter
from pathlib import Path


################################################################################
def process_args():
    docs_args = {
        "main" : "Two modes are available for this standalone calculator: **PocketSphere (PS)** and **Whole (W)**. **PS** calculates the potentials in a spherical volume defined to be a binding pocket, while **W** calculates the potentials of all the volume surrounding the macromolecule.",
        "pdb" : "Path to the PDB structure file of interest.",
        "out" : "Path to the folder where the output potentials should be stored.",
        "apbs" : "Path to the output of APBS for the respective PDB file (this must be done before). An OpenDX file is expected.",
        "rna" : "The target macromolecule of interest is nucleic. If this flag is not provided, the target macromolecule is assumed to be proteic.",
        "whole" : "Activates 'whole mode'. The potentials of all the volume surrounding the macromolecule will be calculated.",
        "default-res" : "Use the default resolution values from settings.py (instead of the default delta values)",
        "radius" : "Radius of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "xcog" : "X coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "ycog" : "Y coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "zcog" : "Z coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.",
    }

    parser = argparse.ArgumentParser(description = docs_args["main"], formatter_class = RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--pdb",    type = Path, required = True, help = docs_args["pdb"])
    parser.add_argument("-o", "--out",    type = Path, required = True, help = docs_args["out"])
    parser.add_argument("-a", "--apbs",   type = Path, default = '',    help = docs_args["apbs"])
    parser.add_argument("-r", "--radius", type = float, default = 1,    help = docs_args["radius"])
    parser.add_argument("-x", "--xcog",   type = float, default = 0,    help = docs_args["xcog"])
    parser.add_argument("-y", "--ycog",   type = float, default = 0,    help = docs_args["ycog"])
    parser.add_argument("-z", "--zcog",   type = float, default = 0,    help = docs_args["zcog"])
    parser.add_argument("-n", "--rna",         action = "store_true",   help = docs_args["rna"])
    parser.add_argument("-w", "--whole",       action = "store_true",   help = docs_args["whole"])
    parser.add_argument("-s", "--default-res", action = "store_true",   help = docs_args["default-res"])
    return parser.parse_args()


################################################################################
