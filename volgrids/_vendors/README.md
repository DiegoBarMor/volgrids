# Vendors

Lightweight packages that are distributed as part of the volgrids package. That way, users don't need to install them separately as a hard dependency.
They aren't kept in the volgrids repo, but instead fetched/packed automatically when releasing to pip.

- If you don't install volgrids via pip, then they will be fetched automatically the first time you run volgrids.
    - Alternatively, run `bash scripts/_prepare.sh` in your local `volgrids` copy to fetch the vendor packages manually.

- Current vendor packages:
    - [FreyaCLI](https://github.com/DiegoBarMor/freyacli) for CLI parsing and management, as well as help string generation and text coloring utilities.
    - [MolUtils](https://github.com/DiegoBarMor/molutils) for utilities dealing with PDB files.
