import os
from pathlib import Path
from abc import ABC, abstractmethod

# //////////////////////////////////////////////////////////////////////////////
class ParamHandler(ABC):
    class InvalidPathError(Exception): pass
    class InvalidParamError(Exception): pass
    class MissingParamError(Exception): pass
    class MissingArgsError(Exception): pass

    # --------------------------------------------------------------------------
    def __init__(self, *params_pos: str, **params_kwd: list[str]):
        self._params_pos = params_pos
        self._params_kwd = params_kwd
        self._help_str: str = ""


    # --------------------------------------------------------------------------
    @abstractmethod
    def assign_globals(self):
        return


    # --------------------------------------------------------------------------
    @property
    @abstractmethod
    def _EXPECTED_CLI_FLAGS() -> dict[str, tuple[str]]:
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @classmethod
    def parse_cli_args(cls, argv) -> tuple[list[str], dict[str, list[str]]]:
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

        for arg in argv:
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
    def _exit_with_help(self, cls_error: type[Exception] = ValueError, err_msg: str = '') -> None:
        if not err_msg: # assume exit code 0 when no error message is provided
            print(f"{self._help_str}")
            exit(0)

        raise cls_error(f"\n{self._help_str}\n\nError while parsing CLI arguments:\n{err_msg}")


    # --------------------------------------------------------------------------
    def _has_params_pos(self) -> bool:
        return len(self._params_pos) > 0


    # --------------------------------------------------------------------------
    def _has_param_kwds(self, *names: str) -> bool:
        return all(name in self._params_kwd for name in names)


    # --------------------------------------------------------------------------
    def _safe_idx(self, lst: list, idx: int, err_msg: str) -> str:
        if idx < 0: idx += len(lst)

        try:
            n_elements = len(lst)
        except TypeError:
            if idx == 0:
                return lst
        else:
            if idx < n_elements:
                return lst[idx]

        self._exit_with_help(err_msg)


    # --------------------------------------------------------------------------
    def _safe_map_value(self, key: str, **values):
        value = values.get(key, None)
        if value is None:
            self._exit_with_help(self.InvalidParamError, f"Invalid parameter '{key}'. Expected one of the following: {', '.join(values.keys())}.")
        return value


    # --------------------------------------------------------------------------
    def _safe_get_param_pos(self, idx: int, err_msg: str = "") -> str:
        return self._safe_idx(self._params_pos, idx, err_msg)


    # --------------------------------------------------------------------------
    def _safe_get_param_kwd_list(self, name: str, min_len: int = 1) -> list[str]:
        if not self._has_param_kwds(name):
            self._exit_with_help(self.MissingParamError, f"The flag '{name}' was not provided.")
        lst = self._params_kwd[name]
        if len(lst) < min_len:
            self._exit_with_help(self.MissingArgsError, f"The flag '{name}' was used but not enough values were provided. At least {min_len} value(s) expected.")
        return lst


    # --------------------------------------------------------------------------
    def _safe_get_param_kwd(self, name: str, required = False, default: str = None) -> str:
        if self._has_param_kwds(name):
            lst = self._params_kwd[name]
            if lst: return lst[0]
            self._exit_with_help(self.MissingArgsError, f"The flag '{name}' was used but no value was provided.")

        if required:
            self._exit_with_help(self.MissingParamError, f"The mandatory flag '{name}' was not provided.")

        return default


    # --------------------------------------------------------------------------
    def _safe_path_file_in(self, path: str) -> Path:
        obj = Path(path)
        if not obj.exists():
            self._exit_with_help(self.InvalidPathError, f"The specified file path '{path}' does not exist.")
        if obj.is_dir():
            self._exit_with_help(self.InvalidPathError, f"The specified file path '{path}' is a folder.")
        return obj


    # --------------------------------------------------------------------------
    def _safe_path_file_out(self, path: str) -> Path:
        obj = Path(path)
        if obj.is_dir():
            self._exit_with_help(self.InvalidPathError, f"The specified file path '{path}' is a folder.")
        os.makedirs(obj.parent, exist_ok = True)
        return obj


    # --------------------------------------------------------------------------
    def _safe_path_folder_out(self, path: str) -> Path:
        obj = Path(path)
        if obj.is_file():
            self._exit_with_help(self.InvalidPathError, f"The specified folder path '{path}' is a file.")
        os.makedirs(obj, exist_ok = True)
        return obj


    # --------------------------------------------------------------------------
    def _safe_kwd_file_in(self, name: str, required = False, default: str = None) -> Path | None:
        path = self._safe_get_param_kwd(name, required, default)
        if path is None: return
        return self._safe_path_file_in(path)


    # --------------------------------------------------------------------------
    def _safe_kwd_file_out(self, name: str, required = False, default: str = None) -> Path | None:
        path = self._safe_get_param_kwd(name, required, default)
        if path is None: return
        return self._safe_path_file_out(path)


    # --------------------------------------------------------------------------
    def _safe_kwd_folder_out(self, name: str, required = False, default: str = None) -> Path | None:
        path = self._safe_get_param_kwd(name, required, default)
        if path is None: return
        return self._safe_path_folder_out(path)


    # --------------------------------------------------------------------------
    def _safe_kwd_float(self, name: str, required = False, default: float = 0.0) -> float:
        val_str: str = self._safe_get_param_kwd(name, required)
        if val_str is None: return default
        try:
            return float(val_str)
        except ValueError:
            self._exit_with_help(self.InvalidParamError, f"The value for the flag '{name}' must be a float. Got '{val_str}' instead.")


# //////////////////////////////////////////////////////////////////////////////
