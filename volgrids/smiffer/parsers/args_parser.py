import os
from pathlib import Path
import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifferArgsParser(vg.ArgsParser):
    def __init__(self):
        super().__init__()
        help_string = '\n'.join((
            "usage: python3 smiffer.py [prot|rna|convert|pack|unpack] [options...]",
            "Available modes:",
            "  prot     - Calculate SMIFs for protein structures.",
            "  rna      - Calculate SMIFs for RNA structures.",
            "  ligand   - Calculate SMIFs for ligand structures. A .chem table must be provided.",
            "Run 'python3 smiffer.py [mode] --help' for more details on each mode.",
            "Running 'python3 smiffer.py' without a valid mode will display this help message.",
        ))

        mode = vg.USER_MODE.lower()
        if mode == "prot":
            sm.CURRENT_MOLTYPE = sm.MolType.PROT
            self._parse_smiffer_calc()
            return

        if mode == "rna":
            sm.CURRENT_MOLTYPE = sm.MolType.RNA
            self._parse_smiffer_calc()
            return

        if mode == "ligand":
            sm.CURRENT_MOLTYPE = sm.MolType.LIGAND
            self._parse_smiffer_calc()
            return

        self.print_exit(-1, help_string)


    # --------------------------------------------------------------------------
    def _parse_smiffer_calc(self) -> None:
        help_string = '\n'.join((
            f"usage: python3 smiffer.py {vg.USER_MODE} [path/input/struct.pdb] [options...]",
            "Available options:",
            "-h, --help                       Show this help message and exit.",
            "-o, --output                     Path to the folder where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --traj                       Path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs                       Path to the output of APBS for the respective structure file (this must be done before). An OpenDX file is expected.",
            "-rxyz", "-ps", "--pocket-sphere  Activate 'pocket sphere' mode by providing the sphere radius and the X, Y, Z coordinates for its center. If not provided, 'whole' mode is assumed.",
            "-b, --table                      Path to a .chem table file to use for ligand mode, or to override the default macromolecules' tables.",
            "-c, --config                     Path to a configuration file with global settings, to override the default settings from config.ini.",
        ))

        fdict = self._get_flags_dict({
            "help"  : ("-h", "--help"),
            "out"   : ("-o", "--output"),
            "traj"  : ("-t", "--traj"),
            "apbs"  : ("-a", "--apbs"),
            "ps"    : ("-rxyz", "-ps", "--pocket-sphere"),
            "table" : ("-b", "--table"),
            "config": ("-c", "--config",)
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_struct = fdict.get(None)
        options_out    = fdict.get("out")
        options_traj   = fdict.get("traj")
        options_apbs   = fdict.get("apbs")
        options_ps     = fdict.get("ps")
        options_table  = fdict.get("table")
        options_config = fdict.get("config")

        if not options_struct:
            self.print_exit(-1, f"{help_string}\nError: No input structure file provided. Provide a path to the structure file as first positional argument.")

        sm.PATH_STRUCTURE  = Path(options_struct[0])

        if not sm.PATH_STRUCTURE.exists():
            self.print_exit(-1, f"{help_string}\nError: The specified structure file '{sm.PATH_STRUCTURE}' does not exist.")

        sm.FOLDER_OUT = Path(options_out[0]) if options_out else sm.PATH_STRUCTURE.parent
        if sm.FOLDER_OUT.is_file():
            self.print_exit(-1, f"{help_string}\nError: The specified output folder '{sm.FOLDER_OUT}' is a file, not a directory.")
        os.makedirs(sm.FOLDER_OUT, exist_ok = True)

        if options_apbs:
            sm.PATH_APBS = Path(options_apbs[0])
            if not sm.PATH_APBS.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified APBS file '{sm.PATH_APBS}' does not exist.")

        if options_ps is not None:
            if len(options_ps) != 4:
                self.print_exit(-1, f"{help_string}\nError: Pocket sphere options must be provided as -rxyz <radius> <x> <y> <z>.")
            try:
                radius = float(options_ps[0])
                x_cog = float(options_ps[1])
                y_cog = float(options_ps[2])
                z_cog = float(options_ps[3])
            except ValueError:
                self.print_exit(-1, f"{help_string}\nError: Pocket sphere options must be numeric values.")
            sm.PS_INFO = (radius, x_cog, y_cog, z_cog)

        if options_traj:
            sm.PATH_TRAJECTORY = Path(options_traj[0])
            if not sm.PATH_TRAJECTORY.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified trajectory file '{sm.PATH_TRAJECTORY}' does not exist.")

        if (sm.CURRENT_MOLTYPE == sm.MolType.LIGAND) and not options_table:
            self.print_exit(-1, f"{help_string}\nError: No table file provided for ligand mode. Use -b or --table to specify the path to the .chem table file.")

        if options_table:
            sm.PATH_TABLE = Path(options_table[0])
            if not sm.PATH_TABLE.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified table file '{sm.PATH_TABLE}' does not exist.")

        if options_config:
            sm.PATH_CONFIG = Path(options_config[0])
            if not sm.PATH_CONFIG.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified config file '{sm.PATH_CONFIG}' does not exist.")


# //////////////////////////////////////////////////////////////////////////////
