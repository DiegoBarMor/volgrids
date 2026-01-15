import sys
import warnings
from pathlib import Path

try:
    import volgrids as vg
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import volgrids as vg

# ------------------------------------------------------------------------------
def help_and_exit(exit_code: int):
    print(
        f"usage: volgrids [smiffer|veins|vgtools] [options...]",
         "Available applications:",
         "  apbs     - Generate raw APBS potential grids for biomolecular structures.",
         "  smiffer  - Calculate SMIFs for biomolecular structures.",
         "  veins    - Calculate VEINS for biomolecular structures.",
         "  vgtools  - Miscellaneous tools for volumetric grids.",
        f"Run 'volgrids [app] --help' for more details on each application.",
        sep = '\n'
    )
    exit(exit_code)

# ------------------------------------------------------------------------------
def run_from_repo():
    vg.PATH_DEFAULT_CONFIG = vg.resolve_path_resource(
        Path(__file__).parent, "config_volgrids.ini"
    )
    main()


# ------------------------------------------------------------------------------
def main():
    argv = sys.argv[1:]
    if not argv: help_and_exit(1)
    if argv[0] in ("-h", "--help"): help_and_exit(0)

    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    app_name = argv[0].lower()
    app_args = argv[1:]

    if app_name == "apbs":
        print(f">>> Launching APBS subprocess for '{app_args[0]}'...", flush = True)
        apbs = vg.APBSSubprocess.run_subprocess(app_args)
        print(f"{apbs.stdout}\n{apbs.stderr}".strip(), flush = True)
        exit(apbs.returncode)

    if app_name == "smiffer":
        import volgrids.smiffer as sm
        sm.AppSmiffer.from_cli(app_args).run()
        exit(0)

    if app_name == "veins":
        import volgrids.veins as ve
        ve.AppVeins.from_cli(app_args).run()
        exit(0)

    if app_name == "vgtools":
        import volgrids.vgtools as vgt
        vgt.AppVGTools.from_cli(app_args).run()
        exit(0)

    help_and_exit(2)


################################################################################
if __name__ == "__main__":
    main()


################################################################################
