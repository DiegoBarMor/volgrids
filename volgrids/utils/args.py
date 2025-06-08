import os, sys, json
from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////
class SmifferArgsParser:
    def __init__(self):
        help_string = '\n'.join((
            "usage: python3 smiffer.py [prot|rna|convert|pack|unpack] [options...]",
            "Available modes:",
            "  prot     - Calculate SMIFs for protein structures.",
            "  rna      - Calculate SMIFs for RNA structures.",
            "  convert  - Convert grid files between formats.",
            "  pack     - Pack multiple grid files into a single CMAP series-file.",
            "  unpack   - Unpack a CMAP series-file into multiple grid files.",
            "  fix-cmap - Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.",
            "Run 'python3 smiffer.py [mode] --help' for more details on each mode.",
            "Running 'python3 smiffer.py' without a valid mode will display this help message.",
        ))

        self.args = sys.argv[1:]
        self.mode = self.next_arg().lower()

        self.path_in: Path = None # "path/input/struct.pdb" SMIFFER
                                  # "path/input/grid.dx"    SMIFF Tools
        self.name: str = None     # name of the input file without extension, e.g. "1abc" for "1abc.pdb"

        ##### SMIFFER
        self.moltype:   str  = None  # protein | nucleic
        self.path_out:  Path = None  # "folder/output/"
        self.path_apbs: Path = None  # "path/input/apbs.pqr.dx"
        self.ps_info: tuple[float, float, float, float] = None
                                     # [radius, x, y, z]
        self.path_traj: Path = None  # "path/input/traj.xtc"
        self.path_meta: Path = None  # automatically set to path_out / f"{self.name}.meta.json"

        self.do_ps:     bool = False # PS mode
        self.do_traj:   bool = False # TRAJ mode


        ##### SMIF Tools
        self.path_dx:   Path = None # "path/input/grid.dx"
        self.path_mrc:  Path = None # "path/output/grid.mrc"
        self.path_cmap: Path = None # "path/output/grid.cmap"


        if self.mode == "prot": # SMIFFER
            self.moltype = "protein"
            self.parse_smiffer_calc()
            return

        if self.mode == "rna": # SMIFFER
            self.moltype = "nucleic"
            self.parse_smiffer_calc()
            return

        if self.mode == "convert": # SMIF Tools
            self.parse_convert()
            return

        if self.mode == "pack": # SMIF Tools
            self.parse_pack()
            return

        if self.mode == "unpack": # SMIF Tools
            self.parse_unpack()
            return

        if self.mode == "fix-cmap": # SMIF Tools
            self.parse_fix_cmap()
            return

        print(help_string)
        exit(-1)


    # --------------------------------------------------------------------------
    def next_arg(self) -> str:
        return self.args.pop(0) if self.args else ''


    # --------------------------------------------------------------------------
    def get_flags_dict(self, flag_names: dict[str, tuple[str]]) -> dict[str, list[str]]:
        """flag_names is a dict where keys are flag identifiers., each associated with a list of aliases for said flag.
        This method then outputs a dict where the keys are the flag identifiers actually found in self.args,
        together with their correspondant values."""
        flags_dict = {None: []} # None is used for options at the start that are not associated with any flag
        alias_to_flagname = {}

        for name,aliases in flag_names.items():
            alias_to_flagname = {**alias_to_flagname, **{alias:name for alias in aliases}}

        current_name = None
        while self.args:
            arg = self.next_arg()
            if arg.lower() in alias_to_flagname: # arg is a flag
                current_name = alias_to_flagname[arg.lower()]
                if current_name in flags_dict:
                    raise ValueError(f"Flag '{current_name}' is defined multiple times in the arguments.")
                flags_dict[current_name] = []

            else: # arg is a flag's option
                flags_dict[current_name].append(arg)

        return flags_dict


    # --------------------------------------------------------------------------
    def parse_smiffer_calc(self) -> None:
        help_string = '\n'.join((
            f"usage: python3 smiffer.py {self.mode} path/input/struct.pdb [options...]",
            "Available options:",
            "-h, --help                       Show this help message and exit.",
            "-o, --output                     Path to the folder where the output SMIFs should be stored. if not provided, the parent folder of the PDB input will be used.",
            "-t, --traj                       Path to a trajectory file (e.g. XTC) supported by MDAnalysis. Activates 'traj' mode: calculate SMIFs for all the frames and save them as a CMAP-series file.",
            "-a, --apbs                       Path to the output of APBS for the respective structure file (this must be done before). An OpenDX file is expected.",
            "-rxyz", "-ps", "--pocket-sphere  Activate 'pocket sphere' mode by providing the sphere radius and the X, Y, Z coordinates for its center. If not provided, 'whole' mode is assumed.",
        ))

        fdict = self.get_flags_dict({
            "help" : ("-h", "--help"),
            "out"  : ("-o", "--output"),
            "traj" : ("-t", "--traj"),
            "apbs" : ("-a", "--apbs"),
            "ps"   : ("-rxyz", "-ps", "--pocket-sphere"),
            # "cav": ("-c", "--cavities"), # [TODO] not implemented yet
        })

        options_in   = fdict.get(None)
        options_out  = fdict.get("out")
        options_traj = fdict.get("traj")
        options_apbs = fdict.get("apbs")
        options_ps   = fdict.get("ps")

        if (fdict.get("help") is not None) or (not options_in):
            if not options_in:
                print(f"{help_string}\nError: No input structure file provided. Provide a path to the structure file as first positional argument.")
                exit(-1)
            print(help_string)
            exit(0)

        self.path_in = Path(options_in[0])
        self.name = self.path_in.stem

        if not self.path_in.exists():
            print(f"{help_string}\nError: The specified structure file '{self.path_in}' does not exist.")
            exit(-1)

        self.path_out = Path(options_out[0]) if options_out \
            else self.path_in.parent
        os.makedirs(self.path_out.parent, exist_ok = True)

        if options_apbs:
            self.path_apbs = Path(options_apbs[0])
            if not self.path_apbs.exists():
                print(f"{help_string}\nError: The specified APBS file '{self.path_apbs}' does not exist.")
                exit(-1)

        if options_ps is not None:
            if len(options_ps) != 4:
                print(f"{help_string}\nError: Pocket sphere options must be provided as -rxyz <radius> <x> <y> <z>.")
                exit(-1)
            try:
                radius = float(options_ps[0])
                x_cog = float(options_ps[1])
                y_cog = float(options_ps[2])
                z_cog = float(options_ps[3])
            except ValueError:
                print(f"{help_string}\nError: Pocket sphere options must be numeric values.")
                exit(-1)

            self.do_ps = True
            self.ps_info = (radius, x_cog, y_cog, z_cog)

        if options_traj:
            self.do_traj = True
            self.path_traj = Path(options_traj[0])
            if not self.path_traj.exists():
                print(f"{help_string}\nError: The specified trajectory file '{self.path_traj}' does not exist.")
                exit(-1)

        self.path_meta = self.path_out / f"{self.name}.meta.json"


    # --------------------------------------------------------------------------
    def parse_convert(self) -> None:
        help_string = '\n'.join((
            "usage: python3 smiffer.py convert path/input/grid.dx [options...]",
            "Available options:",
            "-h, --help  Show this help message and exit.",
            "-d, --dx    Path where to save the converted grid in DX format.",
            "-m, --mrc   Path where to save the converted grid in MRC format.",
            "-c, --cmap  Path where to save the converted grid in CMAP format. The stem of the input file will be used as the CMAP key.",
        ))

        fdict = self.get_flags_dict({
            "help" : ("-h", "--help"),
            "dx"   : ("-d", "--dx"),
            "mrc"  : ("-m", "--mrc"),
            "cmap" : ("-c", "--cmap"),
        })

        options_in   = fdict.get(None)
        options_dx   = fdict.get("dx")
        options_mrc  = fdict.get("mrc")
        options_cmap = fdict.get("cmap")

        if (fdict.get("help") is not None) or (not options_in):
            if not options_in:
                print(f"{help_string}\nError: No input grid file provided. Provide a path to the grid file as first positional argument.")
                exit(-1)
            print(help_string)
            exit(0)

        self.path_in = Path(options_in[0])
        if not self.path_in.exists():
            print(f"{help_string}\nError: The specified grid file '{self.path_in}' does not exist.")
            exit(-1)

        if options_dx:
            self.path_dx = Path(options_dx[0])
            os.makedirs(self.path_dx.parent, exist_ok = True)

        if options_mrc:
            self.path_mrc = Path(options_mrc[0])
            os.makedirs(self.path_mrc.parent, exist_ok = True)

        if options_cmap:
            self.path_cmap = Path(options_cmap[0])
            os.makedirs(self.path_cmap.parent, exist_ok = True)



    # --------------------------------------------------------------------------
    def parse_pack(self) -> None:
        pass


    # --------------------------------------------------------------------------
    def parse_unpack(self) -> None:
        pass


    # --------------------------------------------------------------------------
    def parse_fix_cmap(self) -> None:
        pass


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
            }))


# //////////////////////////////////////////////////////////////////////////////

# import argparse
# def args_tools():
#     docs_args = {
#         "" : "Toolset for post-processing SMIF grids.",
#         "to-mrc" : "Path to the grid file to be converted into MRC format. If input file is a CMAP file, it will be assumed that it contains a single structure, so the first CAMP key is used to access the data.",
#         "to-cmap" : "Path to the grid file to be converted into CMAP format. The stem of the input file will be used as the CMAP key.",
#         "pack" : "Output path of for the packed CMAP series-file, followed by the list of paths to the grid files to be packed.",
#         "unpack" : "Path to the CMAP series-file to be unpacked into several CMAP grid-files.",
#         "fix-cmap" : "Path to the CMAP file to be fixed. This will interpolate all grids into the same resolution.",
#     }

#     parser = argparse.ArgumentParser(description = docs_args[""], formatter_class = argparse.RawDescriptionHelpFormatter)
#     parser.add_argument("-m", "--to-mrc",  type = Path, default = '', help = docs_args["to-mrc"])
#     parser.add_argument("-c", "--to-cmap", type = Path, default = '', help = docs_args["to-cmap"])
#     parser.add_argument("-p", "--pack",    type = Path, default = '', help = docs_args["pack"], nargs = '+')
#     parser.add_argument("-u", "--unpack",  type = Path, default = '', help = docs_args["unpack"], nargs = '+')
#     parser.add_argument("-f", "--fix-cmap",  type = str, default = '', help = docs_args["fix-cmap"], nargs = '+')
#     return parser.parse_args()
