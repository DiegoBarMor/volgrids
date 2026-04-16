import sys
import io
import logging
import warnings
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import MDAnalysis as mda
import pandas as pd

from rnapolis.annotator import  (
    extract_base_interactions,
    handle_input_file,
    read_3d_structure,
    write_csv,
)

# ------------------------------------------------------------------------------
def run_rnapolis(path_pdb: Path, path_csv: Path):
    file = handle_input_file(path_pdb)
    structure3d = read_3d_structure(file, None)
    base_interactions = extract_base_interactions(structure3d)
    structure2d, _ = structure3d.extract_secondary_structure(
        base_interactions, False, False
    )
    write_csv(path_csv, structure2d)


# ------------------------------------------------------------------------------
def get_idxs_canonical(path_csv: Path) -> set[int]:
    df = pd.read_csv(path_csv)
    df = df[
        (df["type"] == "base pair") &
        (df["classification-1"] == "cWW") & (
            (df["classification-2"] == "XIX") |
            (df["classification-2"] == "XX")
        )
    ]
    resid_0 = df["nt1"].apply(lambda x: int(x[3:]))
    resid_1 = df["nt2"].apply(lambda x: int(x[3:]))
    idxs_canonical = set(resid_0) | set(resid_1)
    return idxs_canonical


# ------------------------------------------------------------------------------
def get_mda_universe_quiet(path_pdb) -> mda.Universe:
    buf = io.StringIO()
    logger = logging.getLogger("MDAnalysis")
    old_level = logger.getEffectiveLevel()
    try:
        logger.setLevel(logging.ERROR)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with redirect_stdout(buf), redirect_stderr(buf):
                u = mda.Universe(path_pdb)
    finally:
        logger.setLevel(old_level)
    return u


# ------------------------------------------------------------------------------
def get_idxs_all(path_pdb) -> set[int]:
    u = get_mda_universe_quiet(path_pdb)
    idxs_all = set(int(i) for i in u.residues.resids)
    return idxs_all


# ------------------------------------------------------------------------------
def main():
    run_rnapolis(PATH_PDB, PATH_CSV)
    idxs_canonical = get_idxs_canonical(PATH_CSV)
    idxs_all = get_idxs_all(PATH_PDB)
    idxs_available = idxs_all - idxs_canonical
    print(' '.join(str(i) for i in idxs_available))


################################################################################
if __name__ == "__main__":
    warnings.filterwarnings("ignore", module = "MDAnalysis.*")
    PATH_PDB = Path(sys.argv[1])
    PATH_CSV = Path(sys.argv[2]) if len(sys.argv) > 2 else\
        PATH_PDB.with_suffix(".csv")
    main()


################################################################################
