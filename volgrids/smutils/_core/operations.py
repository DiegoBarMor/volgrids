from pathlib import Path

import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class SMOperations:
    @staticmethod
    def print_resids_nonbp(path_struct: Path) -> None:
        resids_nonbp = su.RNAResids.get_resids_nonbp(path_struct)
        print(resids_nonbp)


# //////////////////////////////////////////////////////////////////////////////
