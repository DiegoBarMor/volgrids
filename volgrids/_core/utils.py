import time
from pathlib import Path

# ------------------------------------------------------------------------------
def resolve_path(path: Path):
    """Resolve the path to the project root directory."""
    project_root = Path(__file__).resolve().parent.parent.parent
    return project_root / path


# //////////////////////////////////////////////////////////////////////////////
class Timer:
    def __init__(self, display = ''):
        if display: print(display, end = ' ', flush = True)
        self.t = None

    def start(self):
        self.t = time.time()

    def end(self):
        if self.t is None:
            raise RuntimeError("Timer not started. Call start() before end().")
        elapsed = time.time() - self.t
        print(f"({int(elapsed // 60)}m {elapsed % 60:.2f}s)", flush = True)


# //////////////////////////////////////////////////////////////////////////////
