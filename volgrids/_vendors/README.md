# Vendors

Lightweight packages can be placed here to be distribitued as part of the `volgrids` package. That way, users don't need to install them separately as a hard dependency.
They aren't kept in the `volgrids` repo, but instead fetched/packed automatically when releasing to pip.
- If you don't install `volgrids` via pip, then run `bash scripts/_prepare_.sh` to fetch the vendors.
- `volgrids` still prioritizes the user installation, if present. Current vendor packages:
    - [FreyaCLI](https://github.com/DiegoBarMor/freyacli) for CLI parsing and management, as well as help string generation and text coloring utilities.
