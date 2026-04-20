from abc import ABC, abstractmethod

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class AppSubcommand(ABC):
    """
    Base class for all subcommands of the volgrids suite of applications.
    They hold a reference to the unique `AppMain` instance, which they can use to access the CLI arguments.
    """

    # --------------------------------------------------------------------------
    def __init__(self, app_main: "vg.AppMain"):
        self.main = app_main


    # --------------------------------------------------------------------------
    @abstractmethod
    def assign_globals(self):
        raise NotImplementedError()


    # --------------------------------------------------------------------------
    @abstractmethod
    def run(self):
        raise NotImplementedError()


# //////////////////////////////////////////////////////////////////////////////
