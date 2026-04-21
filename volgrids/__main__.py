import sys
import warnings
from pathlib import Path

try:
    import volgrids as vg
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import volgrids as vg

# ------------------------------------------------------------------------------
def path_default_config() -> Path | None:
    if not (__package__ == '' or __package__ is None):
        return None
    path_config = Path(__file__).parent.parent / "config_volgrids.ini"
    return path_config if path_config.is_file() else None


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def main():
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    vg.PATH_DEFAULT_CONFIG = path_default_config()
    vg.AppMain(sys.argv).run()


################################################################################
if __name__ == "__main__":
    main()


################################################################################
