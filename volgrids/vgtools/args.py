import os
from pathlib import Path
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class VGToolsArgsParser(vg.ArgsParser):
    def __init__(self):
        super().__init__()
        help_string = '\n'.join((
            "usage: python3 vgtools.py [prot|rna|convert|pack|unpack] [options...]",
            "Available modes:",
            "  convert  - Convert grid files between formats.",
            "  pack     - Pack multiple grid files into a single CMAP series-file.",
            "  unpack   - Unpack a CMAP series-file into multiple grid files.",
            "  fix-cmap - Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.",
            "Run 'python3 vgtools.py [mode] --help' for more details on each mode.",
            "Running 'python3 vgtools.py' without a valid mode will display this help message.",
        ))

        ### Convert
        self.path_dx:   Path = None # "path/input/grid.dx"
        self.path_mrc:  Path = None # "path/output/grid.mrc"
        self.path_cmap: Path = None # "path/output/grid.cmap"

        ### Pack
        self.paths_pack_in: list[Path] = None # list of paths to input grids for packing
        self.path_pack_out: Path = None # "path/output/packed.cmap"

        ### Unpack
        self.path_unpack_in:  Path = None # "path/input/packed.cmap"
        self.path_unpack_out: Path = None # folder where to unpack the grids

        ### Fix CMAP
        self.path_fixcmap_in:  Path = None # "path/input/fix.cmap"
        self.path_fixcmap_out: Path = None # "path/output/fix.cmap"


        if self.mode == "convert":
            self._parse_convert()
            return

        if self.mode == "pack":
            self._parse_pack()
            return

        if self.mode == "unpack":
            self._parse_unpack()
            return

        if self.mode == "fix-cmap":
            self._parse_fix_cmap()
            return

        self.print_exit(-1, help_string)


    # --------------------------------------------------------------------------
    def _parse_convert(self) -> None:
        help_string = '\n'.join((
            "usage: python3 vgtools.py convert path/input/grid.dx [options...]",
            "Available options:",
            "-h, --help  Show this help message and exit.",
            "-d, --dx    Path where to save the converted grid in DX format.",
            "-m, --mrc   Path where to save the converted grid in MRC format.",
            "-c, --cmap  Path where to save the converted grid in CMAP format. The stem of the input file will be used as the CMAP key.",
        ))

        fdict = self._get_flags_dict({
            "help" : ("-h", "--help"),
            "dx"   : ("-d", "--dx"),
            "mrc"  : ("-m", "--mrc"),
            "cmap" : ("-c", "--cmap"),
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_in   = fdict.get(None)
        options_dx   = fdict.get("dx")
        options_mrc  = fdict.get("mrc")
        options_cmap = fdict.get("cmap")

        if not options_in:
            self.print_exit(-1, f"{help_string}\nError: No input grid file provided. Provide a path to the grid file as first positional argument.")

        self.path_in = Path(options_in[0])
        if not self.path_in.exists():
            self.print_exit(-1, f"{help_string}\nError: The specified grid file '{self.path_in}' does not exist.")

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
    def _parse_pack(self) -> None:
        help_string = '\n'.join((
            "usage: python3 vgtools.py pack [options...]",
            "Available options:",
            "-h, --help    Show this help message and exit.",
            "-i, --input   List of paths with the input grids to be packed. At least one grid file must be provided.",
            "-o, --output  Path where to save the packed grid in CMAP format. Must be provided.",
        ))

        fdict = self._get_flags_dict({
            "help"   : ("-h", "--help"),
            "input"  : ("-i", "--input"),
            "output" : ("-o", "--output"),
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_in   = fdict.get("input")
        options_out  = fdict.get("output")

        if not options_in or not options_out:
            self.print_exit(-1, f"{help_string}\nError: Input grids and output path must be provided.")

        self.paths_pack_in = [Path(path) for path in options_in]
        self.path_pack_out = Path(options_out[0])


    # --------------------------------------------------------------------------
    def _parse_unpack(self) -> None:
        help_string = '\n'.join((
            "usage: python3 vgtools.py unpack [options...]",
            "Available options:",
            "-h, --help    Show this help message and exit.",
            "-i, --input   Path to the CMAP series-file to be unpacked. Must be provided.",
            "-o, --output  Path where to save the unpacked grids. If not provided, the parent folder of the input packed file will be used.",
        ))

        fdict = self._get_flags_dict({
            "help"   : ("-h", "--help"),
            "input"  : ("-i", "--input"),
            "output" : ("-o", "--output"),
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_in   = fdict.get("input")
        options_out  = fdict.get("output")

        if not options_in:
            self.print_exit(-1, f"{help_string}\nError: Input packed file must be provided.")

        self.path_unpack_in = Path(options_in[0])
        if not self.path_unpack_in.exists():
            self.print_exit(-1, f"{help_string}\nError: The specified packed file '{self.path_unpack_in}' does not exist.")

        self.path_unpack_out = Path(options_out[0]) if options_out else self.path_unpack_in.parent
        if self.path_unpack_out.is_file():
            self.print_exit(-1, f"{help_string}\nError: The specified output folder '{self.path_unpack_out}' is a file, not a directory.")


    # --------------------------------------------------------------------------
    def _parse_fix_cmap(self) -> None:
        help_string = '\n'.join((
            "usage: python3 vgtools.py fix-cmap [options...]",
            "Available options:",
            "-h, --help    Show this help message and exit.",
            "-i, --input   Path to the CMAP file to be fixed. Must be provided.",
            "-o, --output  Path where to save the fixed CMAP file. Must be provided.",
        ))

        fdict = self._get_flags_dict({
            "help"   : ("-h", "--help"),
            "input"  : ("-i", "--input"),
            "output" : ("-o", "--output"),
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_in   = fdict.get("input")
        options_out  = fdict.get("output")

        if not options_in:
            self.print_exit(-1, f"{help_string}\nError: No input CMAP file provided.")

        if not options_out:
            self.print_exit(-1, f"{help_string}\nError: No output path provided for the fixed CMAP file.")

        self.path_fixcmap_in  = Path(options_in[0])
        if not self.path_fixcmap_in.exists():
            self.print_exit(-1, f"{help_string}\nError: The specified CMAP file '{self.path_fixcmap_in}' does not exist.")

        self.path_fixcmap_out = Path(options_out[0])
        os.makedirs(self.path_fixcmap_out.parent, exist_ok = True)


# //////////////////////////////////////////////////////////////////////////////
