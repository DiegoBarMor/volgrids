# Vendors

Lightweight packages can be placed here to be distribitued as part of the volgrids package. That way, users don't need to install them separately as a hard dependency.
They aren't kept in the volgrids repo, but instead fetched/packed automatically when releasing to pip.

- If you don't install volgrids via pip, then they will be fetched automatically the first time you run volgrids.
    - Alternatively, run `bash volgrids/_prepare.sh` in your local `volgrids` copy to fetch the vendor packages manually.
    - volgrids still prioritizes the vendors' pip installation, if present.

- Current vendor packages:
    - [FreyaCLI](https://github.com/DiegoBarMor/freyacli) for CLI parsing and management, as well as help string generation and text coloring utilities.
