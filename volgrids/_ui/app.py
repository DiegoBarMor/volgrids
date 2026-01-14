from abc import ABC, abstractmethod
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class App(ABC):
    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        handler = self._CLASS_PARAM_HANDLER(*args, **kwargs)
        handler.assign_globals()
        self.load_configs()


    # --------------------------------------------------------------------------
    @classmethod
    def from_cli(cls, argv):
        params_pos, params_kwd = cls._CLASS_PARAM_HANDLER.parse_cli_args(argv)
        return cls(*params_pos, **params_kwd)


    # --------------------------------------------------------------------------
    @classmethod
    def from_str_argv(cls, str_argv: str):
        argv = str_argv.split()
        params_pos, params_kwd = cls._CLASS_PARAM_HANDLER.parse_cli_args(argv)
        return cls(*params_pos, **params_kwd)


    # --------------------------------------------------------------------------
    def load_configs(self):
        if vg.PATH_DEFAULT_CONFIG is not None:
            self._load_config_file(vg.PATH_DEFAULT_CONFIG)

        if vg.PATH_CUSTOM_CONFIG is not None:
            self._load_config_file(vg.PATH_CUSTOM_CONFIG)


    # --------------------------------------------------------------------------
    def _load_config_file(self, path_config: Path) -> None:
        parser = vg.ParserConfig(path_config)
        scope_dependencies = self._import_config_dependencies()

        for section, scope_module in self.CONFIG_MODULES.items():
            parser.apply_config(
                section = section,
                scope_module = scope_module.__dict__,
                scope_dependencies = scope_dependencies,
                valid_config_keys = scope_module.__config_keys__.copy(),
            )


    # --------------------------------------------------------------------------
    @property
    @abstractmethod
    def CONFIG_MODULES() -> dict[str, dict[str, any]]:
        """class variable that links the pertinent keys from the config file, to the module to be configured.
        For example: `{"VOLGRIDS": vg, "SMIFFER": sm}`, assuming `import volgrids.vgrids as vg` and `import volgrids.smiffer as sm`."""
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @property
    @abstractmethod
    def _CLASS_PARAM_HANDLER() -> type["vg.ParamHandler"]:
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @abstractmethod
    def run(self):
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @abstractmethod
    def _import_config_dependencies(self) -> dict[str, any]:
        """Import any modules that can be used by the config file.
        For example, if a key-value pair is VALUE=np.inf, you should import numpy as np inside this method.
        Return the scope with the imported modules, for example by adding "return locals()" at the end.
        """
        raise NotImplementedError()


# //////////////////////////////////////////////////////////////////////////////
