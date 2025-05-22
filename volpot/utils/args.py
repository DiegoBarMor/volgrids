import argparse
from pathlib import Path

################################################################################
def args_smiffer():
    docs_args = {
        "" : "Two modes are available for this standalone calculator: **PocketSphere (PS)** and **Whole (W)**. **PS** calculates the potentials in a spherical volume defined to be a binding pocket, while **W** calculates the potentials of all the volume surrounding the macromolecule.",
        "pdb" : "Path to the PDB structure file of interest.",
        "out" : "Path to the folder where the output potentials should be stored.",
        "traj" : "Path to the trajectory file for TRAJ mode.",
        "apbs" : "Path to the output of APBS for the respective PDB file (this must be done before). An OpenDX file is expected.",
        "radius" : "Radius of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "xcog" : "X coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "ycog" : "Y coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "zcog" : "Z coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.",
        "rna" : "The target macromolecule of interest is nucleic. If this flag is not provided, the target macromolecule is assumed to be proteic.",
        "whole" : "Activates 'whole mode'. The potentials of all the volume surrounding the macromolecule will be calculated.",
        "default-res" : "Use the default resolution values from settings.py (instead of the default delta values)",
    }

    parser = argparse.ArgumentParser(description = docs_args[""], formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--pdb",    type = Path, required = True, help = docs_args["pdb"])
    parser.add_argument("-o", "--out",    type = Path, required = True, help = docs_args["out"])
    parser.add_argument("-t", "--traj",   type = Path, default = '',    help = docs_args["traj"])
    parser.add_argument("-a", "--apbs",   type = Path, default = '',    help = docs_args["apbs"])
    parser.add_argument("-r", "--radius", type = float, default = 1,    help = docs_args["radius"])
    parser.add_argument("-x", "--xcog",   type = float, default = 0,    help = docs_args["xcog"])
    parser.add_argument("-y", "--ycog",   type = float, default = 0,    help = docs_args["ycog"])
    parser.add_argument("-z", "--zcog",   type = float, default = 0,    help = docs_args["zcog"])
    parser.add_argument("-n", "--rna",         action = "store_true",   help = docs_args["rna"])
    parser.add_argument("-w", "--whole",       action = "store_true",   help = docs_args["whole"])
    parser.add_argument("-s", "--default-res", action = "store_true",   help = docs_args["default-res"])
    return parser.parse_args()


def args_tools():
    docs_args = {
        "" : "Toolset for post-processing SMIF grids.",
        "to-mrc" : "Path to the grid file to be converted into MRC format. If input file is a CMAP file, it will be assumed that it contains a single structure, so the first CAMP key is used to access the data.",
        "to-cmap" : "Path to the grid file to be converted into CMAP format. The stem of the input file will be used as the CMAP key.",
        "pack" : "Output path of for the packed CMAP series-file, followed by the list of paths to the grid files to be packed.",
        "unpack" : "Path to the CMAP series-file to be unpacked into several CMAP grid-files.",
    }

    parser = argparse.ArgumentParser(description = docs_args[""], formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--to-mrc",  type = Path, default = '', help = docs_args["to-mrc"])
    parser.add_argument("-c", "--to-cmap", type = Path, default = '', help = docs_args["to-cmap"])
    parser.add_argument("-p", "--pack",    type = Path, default = '', help = docs_args["pack"], nargs = '+')
    parser.add_argument("-u", "--unpack",  type = Path, default = '', help = docs_args["unpack"])
    return parser.parse_args()


################################################################################
