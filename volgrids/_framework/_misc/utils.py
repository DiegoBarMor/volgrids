from pathlib import Path

# ------------------------------------------------------------------------------
def resolve_path(path: Path):
    """Resolve the path to the project root directory."""
    project_root = Path(__file__).resolve().parent.parent.parent
    return project_root / path


# ------------------------------------------------------------------------------
