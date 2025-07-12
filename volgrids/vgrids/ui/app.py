from abc import ABC, abstractmethod

import volgrids.vgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class App(ABC):
    def __init__(self, *args, **kwargs):
        handler = self._CLASS_PARAM_HANDLER(*args, **kwargs)
        handler.assign_globals()


    # --------------------------------------------------------------------------
    @classmethod
    def from_cli(cls):
        params_pos, params_kwd = cls._CLASS_PARAM_HANDLER.parse_cli_args()
        return cls(*params_pos, **params_kwd)


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
