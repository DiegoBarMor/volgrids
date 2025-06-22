import os, json
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

        self.moltype: sm.MolType = sm.MolType.NONE # type of the molecule, e.g. PROT, RNA, LIGAND
        self.path_in: Path = None # "path/input/struct.pdb"
        self.name: str = None     # name of the input file without extension, e.g. "1abc" for "1abc.pdb"

        self.ps_info: tuple[float, float, float, float] = None # [radius, x, y, z]
        self.path_out:   Path = None # "folder/output/"
        self.path_apbs:  Path = None # "path/input/apbs.pqr.dx"
        self.path_traj:  Path = None # "path/input/traj.xtc"
        self.path_table: Path = None # "path/input/table.chem"
        self.path_meta:  Path = None # automatically set to path_out / f"{self.name}.meta.json"

        self.do_ps:   bool = False # PS mode
        self.do_traj: bool = False # TRAJ mode


        if self.mode == "prot":
            self.moltype = sm.MolType.PROT
            self._parse_smiffer_calc()
            return

        if self.mode == "rna":
            self.moltype = sm.MolType.RNA
            self._parse_smiffer_calc()
            return

        if self.mode == "ligand":
            self.moltype = sm.MolType.LIGAND
            self._parse_smiffer_calc()
            return

        self.print_exit(-1, help_string)


    # --------------------------------------------------------------------------
    def save_metadata(self):
        if self.path_meta is None:
            print("No metadata path specified. Metadata will not be saved.")
            return

        with open(self.path_meta, 'w') as file:
            file.write(json.dumps({
                "mode": self.mode,
                "input":  str(self.path_in),
                "output": str(self.path_out),
                "apbs":   str(self.path_apbs),
                "pocket_sphere": self.ps_info,
                "trajectory": str(self.path_traj),
                "table": str(self.path_table),
            }))


    # --------------------------------------------------------------------------
    def _parse_smiffer_calc(self) -> None:
        help_string = '\n'.join((
            f"usage: python3 smiffer.py {self.mode} path/input/struct.pdb [options...]",
            "Available options:",
            "-h, --help                       Show this help message and exit.",
            "-o, --output                     Path to the folder where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --traj                       Path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs                       Path to the output of APBS for the respective structure file (this must be done before). An OpenDX file is expected.",
            "-rxyz", "-ps", "--pocket-sphere  Activate 'pocket sphere' mode by providing the sphere radius and the X, Y, Z coordinates for its center. If not provided, 'whole' mode is assumed.",
            "-b, --table                      Path to a .chem table file to use for ligand mode, or to override the default macromolecules' tables.",
        ))

        fdict = self._get_flags_dict({
            "help" : ("-h", "--help"),
            "out"  : ("-o", "--output"),
            "traj" : ("-t", "--traj"),
            "apbs" : ("-a", "--apbs"),
            "ps"   : ("-rxyz", "-ps", "--pocket-sphere"),
            "table": ("-b", "--table"),
            # "cav": ("-c", "--cavities"), # [TODO] not implemented yet
            "debug": ("--debug",)
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_in    = fdict.get(None)
        options_out   = fdict.get("out")
        options_traj  = fdict.get("traj")
        options_apbs  = fdict.get("apbs")
        options_ps    = fdict.get("ps")
        options_table = fdict.get("table")

        if not options_in:
            self.print_exit(-1, f"{help_string}\nError: No input structure file provided. Provide a path to the structure file as first positional argument.")

        self.path_in = Path(options_in[0])
        self.name = self.path_in.stem

        if not self.path_in.exists():
            self.print_exit(-1, f"{help_string}\nError: The specified structure file '{self.path_in}' does not exist.")

        self.path_out = Path(options_out[0]) if options_out else self.path_in.parent
        if self.path_out.is_file():
            self.print_exit(-1, f"{help_string}\nError: The specified output folder '{self.path_out}' is a file, not a directory.")

        os.makedirs(self.path_out, exist_ok = True)

        if options_apbs:
            self.path_apbs = Path(options_apbs[0])
            if not self.path_apbs.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified APBS file '{self.path_apbs}' does not exist.")

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

            self.do_ps = True
            self.ps_info = (radius, x_cog, y_cog, z_cog)

        if options_traj:
            self.do_traj = True
            self.path_traj = Path(options_traj[0])
            if not self.path_traj.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified trajectory file '{self.path_traj}' does not exist.")

        if (self.moltype == sm.MolType.LIGAND) and not options_table:
            self.print_exit(-1, f"{help_string}\nError: No table file provided for ligand mode. Use -b or --table to specify the path to the .chem table file.")

        if options_table:
            self.path_table = Path(options_table[0])
            if not self.path_table.exists():
                self.print_exit(-1, f"{help_string}\nError: The specified table file '{self.path_table}' does not exist.")

        self.path_meta = self.path_out / f"{self.name}.meta.json"

        if fdict.get("debug"):
            self._get_debug_vars(fdict["debug"])


# //////////////////////////////////////////////////////////////////////////////
