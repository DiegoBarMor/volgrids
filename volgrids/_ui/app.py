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
    def load_configs(self) -> None:
        self._load_config(vg.PATH_DEFAULT_CONFIG)
        for path_config in vg.PATHS_CUSTOM_CONFIG:
            self._load_config(path_config)
        self._load_config(vg.STR_CUSTOM_CONFIG, is_file = False)


    # --------------------------------------------------------------------------
    def _load_config(self, config: Path | str, is_file = True) -> None:
        if is_file:
            if config is None: return
            parser = vg.ParserConfig.from_file(config)
        else:
            if not config.strip(): return
            parser = vg.ParserConfig(config)

        all_known_keys = set(k for scope_module in self.CONFIG_MODULES for k in scope_module.__config_keys__)
        for scope_module in self.CONFIG_MODULES:
            parser.apply_config(
                scope_module = scope_module.__dict__,
                this_module_keys = scope_module.__config_keys__,
                all_known_keys = all_known_keys
            )


    # --------------------------------------------------------------------------
    @property
    @abstractmethod
    def CONFIG_MODULES() -> tuple:
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


# //////////////////////////////////////////////////////////////////////////////
