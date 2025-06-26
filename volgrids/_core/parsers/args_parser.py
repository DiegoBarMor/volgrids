import sys

# //////////////////////////////////////////////////////////////////////////////
class ArgsParser:
    """A class to parse command line arguments for the VolGrids-based applications.
    It assumes that the application is run with the following structure:\n
    `python3 app.py [mode] [options...]`\n
    where `[mode]` is one of multiple available modes
    and `[options...]` are flags/options for the respective mode."""

    # --------------------------------------------------------------------------
    def __init__(self):
        self.args = sys.argv[1:]
        self.mode = self._next_arg().lower()


    # --------------------------------------------------------------------------
    def print_exit(self, exit_code: int, message: str):
        print(message)
        exit(exit_code)


    # --------------------------------------------------------------------------
    def _next_arg(self) -> str:
        return self.args.pop(0) if self.args else ''


    # --------------------------------------------------------------------------
    def _get_flags_dict(self, flag_names: dict[str, tuple[str]]) -> dict[str, list[str]]:
        """flag_names is a dict where keys are flag identifiers., each associated with a list of aliases for said flag.
        This method then outputs a dict where the keys are the flag identifiers actually found in self.args,
        together with their correspondant values."""
        flags_dict = {None: []} # None is used for options at the start that are not associated with any flag
        alias_to_flagname = {}

        for name,aliases in flag_names.items():
            alias_to_flagname = {**alias_to_flagname, **{alias:name for alias in aliases}}

        current_name = None
        while self.args:
            arg = self._next_arg()
            if arg.lower() in alias_to_flagname: # arg is a flag
                current_name = alias_to_flagname[arg.lower()]
                if current_name in flags_dict:
                    raise ValueError(f"Flag '{current_name}' is defined multiple times in the arguments.")
                flags_dict[current_name] = []

            else: # arg is a flag's option
                flags_dict[current_name].append(arg)

        return flags_dict


# //////////////////////////////////////////////////////////////////////////////
