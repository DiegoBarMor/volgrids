import os
from pathlib import Path
import volgrids as vg
import volgrids.veins as ve

# //////////////////////////////////////////////////////////////////////////////
class VeinsArgsParser(vg.ArgsParser):
    def __init__(self):
        super().__init__()
        help_string = '\n'.join((
            "usage: python3 veins.py [mode] [options...]",
            "Available modes:",
            "  energies - Generate grids to visually represent in space the interaction energies of a molecular system.",
            # "  force  - ", # TODO: Implement force mode
            "Run 'python3 veins.py [mode] --help' for more details on each mode.",
            "Running 'python3 veins.py' without a valid mode will display this help message.",
        ))

        mode = vg.USER_MODE.lower()
        if mode == "energies":
            self._parse_energies()
            return

        # if mode == "force":
        #     self._parse_forces()
        #     return

        self.print_exit(-1, help_string)


    # --------------------------------------------------------------------------
    def _parse_energies(self) -> None:
        help_string = '\n'.join((
            "usage: python3 veins.py energy [path/input/structure.pdb] [path/input/energies.csv] [options...]",
            "Available options:",
            "-h, --help       Show this help message and exit.",
            "-o, --output     Path to the folder where the output SMIFs should be stored. If not provided, the parent folder of the input structure file will be used.",
            "-t, --trajectory Path to a trajectory file (e.g. XTC) supported by MDAnalysis. In this case, the energies CSV file contains an energy column for each frame. The header of such columns must start with 'frame'.",
            "-c, --cutoff     Energies below this cutoff will be ignored. Default value: 1e-3.",
        ))

        fdict = self._get_flags_dict({
            "help" : ("-h", "--help"),
            "out"  : ("-o", "--output"),
            "traj" : ("-t", "--trajectory"),
            "cut"  : ("-c", "--cutoff"),
        })

        if fdict.get("help") is not None:
            self.print_exit(0, help_string)

        options_in   = fdict.get(None)
        options_out  = fdict.get("out")
        options_cut  = fdict.get("cut")
        options_traj = fdict.get("traj")

        if len(options_in) < 2:
            self.print_exit(-1, f"{help_string}\nError: You must provide both the input structure file and the energies CSV file.")

        ve.PATH_STRUCTURE = Path(options_in[0])
        if not ve.PATH_STRUCTURE.exists():
            self.print_exit(-1, f"Error: The input structure file '{ve.PATH_STRUCTURE}' does not exist.")

        ve.PATH_ENERGIES_CSV = Path(options_in[1])
        if not ve.PATH_ENERGIES_CSV.exists():
            self.print_exit(-1, f"Error: The energies CSV file '{ve.PATH_ENERGIES_CSV}' does not exist.")

        ve.FOLDER_OUT = Path(options_out[0]) if options_out else ve.PATH_STRUCTURE.parent
        if ve.FOLDER_OUT.is_file():
            self.print_exit(-1, f"{help_string}\nError: The specified output folder '{ve.FOLDER_OUT}' is a file, not a directory.")
        os.makedirs(ve.FOLDER_OUT, exist_ok = True)

        if options_traj:
            ve.PATH_TRAJECTORY = Path(options_traj[0])
            if not ve.PATH_TRAJECTORY.exists():
                self.print_exit(-1, f"Error: The trajectory file '{ve.PATH_TRAJECTORY}' does not exist.")

        if options_cut:
            try:
                ve.ENERGY_CUTOFF = float(options_cut[0])
            except ValueError:
                self.print_exit(-1, f"{help_string}\nError: The cutoff value must be a valid float number. Provided: '{options_cut[0]}'")


# //////////////////////////////////////////////////////////////////////////////
