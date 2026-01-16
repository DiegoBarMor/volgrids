from pathlib import Path

# ------------------------------------------------------------------------------
def resolve_path_package(path_to_resolve: str | Path):
    """Resolve the path to the package root directory."""
    root_package = Path(__file__).resolve().parent.parent
    return root_package / path_to_resolve


# ------------------------------------------------------------------------------
