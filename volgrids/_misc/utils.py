import subprocess
from pathlib import Path

# ------------------------------------------------------------------------------
def resolve_path_package(path_to_resolve: str | Path):
    """Resolve the path to the package root directory."""
    root_package = Path(__file__).resolve().parent.parent
    return root_package / path_to_resolve


# ------------------------------------------------------------------------------
def assert_vendors():
    """
    Automatic fetching of the vendors first time volgrids is executed,
    in case volgrids isn't downloaded via pip.
    """
    dir_vendors = resolve_path_package("_vendors")
    names_vendors = ("freyacli", "molutils")
    if all((dir_vendors / name).is_dir() for name in names_vendors):
        return

    print(">>> Fetching vendor dependencies for the first time...")

    path_sh = resolve_path_package("_vendors/fetch_vendors.sh")
    proc = subprocess.run(
        ["/bin/bash", str(path_sh)],
        capture_output = True, text = True
    )
    if proc.returncode == 0: return

    raise RuntimeError('\n'.join((
        f"Failed to fetch vendors (code={proc.returncode}):",
        proc.stderr or "<empty_stderr>",
        proc.stdout or "<empty_stdout>",
    )))


# ------------------------------------------------------------------------------
