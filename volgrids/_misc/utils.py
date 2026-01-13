from pathlib import Path

# ------------------------------------------------------------------------------
def resolve_path_package(path_to_resolve: str | Path):
    """Resolve the path to the package root directory."""
    root_package = Path(__file__).resolve().parent.parent
    return root_package / path_to_resolve


# ------------------------------------------------------------------------------
def resolve_path_resource(path_py: str, path_to_resolve: str | Path):
    """Resolve the path to a resource, given the path to the calling python file."""
    parent_py = Path(path_py).resolve().parent
    return parent_py / path_to_resolve


# ------------------------------------------------------------------------------
