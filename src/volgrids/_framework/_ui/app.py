from abc import ABC, abstractmethod
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class App(ABC):
    PATH_DEFAULT_CONFIG = vg.resolve_path("config.ini")

    # --------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        handler = self._CLASS_PARAM_HANDLER(*args, **kwargs)
        handler.assign_globals()
        self.load_configs(vg.PATH_CUSTOM_CONFIG)


    # --------------------------------------------------------------------------
    @classmethod
    def from_cli(cls):
        params_pos, params_kwd = cls._CLASS_PARAM_HANDLER.parse_cli_args()
        return cls(*params_pos, **params_kwd)


    # --------------------------------------------------------------------------
    def load_configs(self, path_custom: Path | None):
        self._load_config_file(self.PATH_DEFAULT_CONFIG, is_default = True)
        if path_custom is None: return
        self._load_config_file(path_custom, is_default = False)


    # --------------------------------------------------------------------------
    def _load_config_file(self, path_config: Path, is_default: bool) -> None:
        parser = vg.ParserConfig(path_config)
        scope_dependencies = self._import_config_dependencies()

        for section, scope_module in self.CONFIG_MODULES.items():
            parser.apply_config(
                section = section,
                scope_module = scope_module.__dict__,
                scope_dependencies = scope_dependencies,
                valid_config_keys = scope_module.__config_keys__.copy(),
                all_configs_mandatory = is_default
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
