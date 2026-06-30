from pathlib import Path

import volgrids as vg
import volgrids.smiffer as smf
import volgrids.smutils as sut
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppBoxes(vg.AppSubcommand):
    # -------------------------------------------------------------------------- UI SECTION
    def run(self):
        operation = self.main.subcommands.pop(0)
        if operation == "info": return self.app_run_info()
        if operation == "size": return self.app_run_size()
        raise ValueError(f"Unknown operation: {operation}")


    # --------------------------------------------------------------------------
    def app_run_info(self):
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)

        vals = sut.AppBoxes.info(path_in)
        print(*(f"{v:.3f}" for v in vals))


    # --------------------------------------------------------------------------
    def app_run_size(self):
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)

        vals = sut.AppBoxes.info(path_in)
        print(*(
            f"{v1-v0:.3f}" for v0,v1 in zip(vals[::2], vals[1::2])
        ))


    # -------------------------------------------------------------------------- LOGIC SECTION
    @staticmethod
    def info(path_pdb: Path) -> tuple[float]:
        return vg.BoxInfo.from_box(
            smf.MoleculeManager(path_pdb).box
        ).values()


# //////////////////////////////////////////////////////////////////////////////
