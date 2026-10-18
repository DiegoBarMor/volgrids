import os
import tempfile
from typing import IO


def handle_input_file(path) -> IO[str]:
    _, ext = os.path.splitext(path)
    file = tempfile.NamedTemporaryFile("wt+", suffix=ext)
    with open(path) as f:
        file.write(f.read())
        file.seek(0)
    return file
