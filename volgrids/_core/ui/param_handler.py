import os, sys
from pathlib import Path
from abc import ABC, abstractmethod

# //////////////////////////////////////////////////////////////////////////////
class ParamHandler(ABC):
    _EXPECTED_CLI_FLAGS: dict[str, tuple[str]] = None # defined in subclasses

    def __init__(self, *params_pos: str, **params_kwd: list[str]):
        if not (params_pos or params_kwd):
            raise ValueError("No parameters have been provided to the ParamHandler.")

        self._params_pos = params_pos
        self._params_kwd = params_kwd
        self._help_str: str = ""


    # --------------------------------------------------------------------------
    @abstractmethod
    def assign_globals(self):
        return


    # --------------------------------------------------------------------------
    @classmethod
    def parse_cli_args(cls) -> tuple[list[str], dict[str, list[str]]]:
        """_EXPECTED_CLI_FLAGS is a dict where keys are flag identifiers, each associated with a list of aliases for said flag.
        This method then outputs a dict where the keys are the flag identifiers actually found in self._args,
        together with their correspondant values."""

        if cls._EXPECTED_CLI_FLAGS is None:
            raise NotImplementedError("The _EXPECTED_CLI_FLAGS attribute must be defined in the subclass.")

        alias_to_flagname = {}
        for name,aliases in cls._EXPECTED_CLI_FLAGS.items():
            alias_to_flagname = {**alias_to_flagname, **{alias:name for alias in aliases}}

        current_name = '' # '' is used for options at the start that are not associated with any flag
        params_kwd: dict[str, list[str]] = {current_name: []}
        cli_args = sys.argv[1:]

        while cli_args:
            arg = cli_args.pop(0)
            if arg.lower() in alias_to_flagname: # arg is a flag
                current_name = alias_to_flagname[arg.lower()]
                if current_name in params_kwd:
                    raise ValueError(f"Flag '{current_name}' is defined multiple times in the arguments.")
                params_kwd[current_name] = []

            else: # arg is a flag's option
                params_kwd[current_name].append(arg)

        params_pos = params_kwd.pop('')
        return params_pos, params_kwd


    # --------------------------------------------------------------------------
    def _set_help_str(self, *lines: str) -> None:
        self._help_str = '\n'.join(lines)


    # --------------------------------------------------------------------------
    def _exit_with_help(self, exit_code: int, err_msg: str = '') -> None:
        print(f"{self._help_str}\n\nError: {err_msg}")
        exit(exit_code)


    # --------------------------------------------------------------------------
    def _has_params_pos(self) -> bool:
        return len(self._params_pos) > 0


    # --------------------------------------------------------------------------
    def _has_param_kwds(self, *names: str) -> list[str]:
        return all(name in self._params_kwd for name in names)


    # --------------------------------------------------------------------------
    def _safe_idx(self, lst: list, idx: int, err_msg: str) -> str:
        if idx < 0: idx += len(lst)
        if len(lst) > idx:
            return lst[idx]
        self._exit_with_help(-1, err_msg)


    # --------------------------------------------------------------------------
    def _safe_map_value(self, key: str, **values):
        value = values.get(key.lower(), None)
        if value is None:
            self._exit_with_help(-1, f"Invalid parameter '{key}'. Expected one of the following: {', '.join(values.keys())}.")
        return value


    # --------------------------------------------------------------------------
    def _safe_get_param_pos(self, idx: int, err_msg: str = "") -> str:
        return self._safe_idx(self._params_pos, idx, err_msg)


    # --------------------------------------------------------------------------
    def _safe_get_param_kwd_list(self, name: str, err_msg: str = None) -> list[str]:
        if not self._has_param_kwds(name):
            self._exit_with_help(-1, err_msg if err_msg else f"The flag '{name}' was not provided.")
        return self._params_kwd[name]


    # --------------------------------------------------------------------------
    def _safe_get_param_kwd(self, name: str, idx: int, err_msg: str = None) -> str:
        lst = self._safe_get_param_kwd_list(name, err_msg)
        return self._safe_idx(lst, idx, err_msg if err_msg else f"The flag '{name}' was used but not enough values were provided.")


    # --------------------------------------------------------------------------
    def _safe_path_file_in(self, path: str) -> Path:
        obj = Path(path)
        if not obj.exists():
            self._exit_with_help(-1, f"The specified file path '{path}' does not exist.")
        if obj.is_dir():
            self._exit_with_help(-1, f"The specified file path '{path}' is a folder.")
        return obj


    # --------------------------------------------------------------------------
    def _safe_path_file_out(self, path: str) -> Path:
        obj = Path(path)
        if obj.is_dir():
            self._exit_with_help(-1, f"The specified file path '{path}' is a folder.")
        os.makedirs(obj.parent, exist_ok = True)
        return obj


    # --------------------------------------------------------------------------
    def _safe_path_folder_out(self, path: str) -> Path:
        obj = Path(path)
        if obj.is_file():
            self._exit_with_help(-1, f"The specified folder path '{path}' is a file.")
        os.makedirs(obj, exist_ok = True)
        return obj


    # --------------------------------------------------------------------------
    def _safe_kwd_file_in(self, name: str) -> Path:
        return self._safe_path_file_in(
            self._safe_get_param_kwd(name, 0)
        )


    # --------------------------------------------------------------------------
    def _safe_kwd_file_out(self, name: str) -> Path:
        return self._safe_path_file_out(
            self._safe_get_param_kwd(name, 0)
        )


# //////////////////////////////////////////////////////////////////////////////
