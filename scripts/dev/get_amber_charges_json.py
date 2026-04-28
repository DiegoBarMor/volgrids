import json
import pandas as pd
from pathlib import Path

NBASE_PUR = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"]
NBASE_PYR = ["N1", "C2", "N3", "C4", "C5", "C6"]

RNA_PHOSP = ["P", "OP1", "OP2", "O5'", "O3'"]
RNA_SUGAR = ["C1'", "C2'", "C3'", "C4'", "O4'", "O2'", "HO2'", "H1'", "H2'", "H3'", "H4'"]
RNA_NBASE = {
    "U": NBASE_PYR + ["O2", "O4", "H3", "H5", "H6"],
    "C": NBASE_PYR + ["O2", "N4", "H41", "H42", "H5", "H6"],
    "A": NBASE_PUR + ["N6", "H61", "H62", "H2", "H8"],
    "G": NBASE_PUR + ["O6", "N2", "H21", "H22", "H1", "H8"],
}

DNA_PHOSP = RNA_PHOSP.copy()
DNA_SUGAR = ["C1'", "C2'", "C3'", "C4'", "O4'", "H1'", "H2'", "H2''", "H3'", "H4'"]
DNA_NBASE = {
    "DT": NBASE_PYR + ["O2", "O4", "C7", "H71", "H72", "H73", "H3", "H6"],
    "DC": RNA_NBASE["C"].copy(),
    "DA": RNA_NBASE["A"].copy(),
    "DG": RNA_NBASE["G"].copy(),
}

# ------------------------------------------------------------------------------
def csv_to_json(path_csv: Path, path_json: Path):
    df = pd.read_csv(path_csv)
    resnames = df["resname"].unique()
    out = {}
    for resname in resnames:
        df_res = df[df["resname"] == resname]
        out[resname] = dict(zip(
            df_res["atomname"],
            list(zip(df_res["charge"], df_res["radius"]))
        ))

    with open(path_json, 'w') as file:
        json.dump(out, file, indent = 4)


# ------------------------------------------------------------------------------
def save_coarse_charges(path_json_in: Path, path_json_out: Path, is_rna: bool):
    if is_rna:
        atoms_phosp = RNA_PHOSP
        atoms_sugar = RNA_SUGAR
        atoms_nbase = RNA_NBASE
        resnames = ("U", "C", "A", "G")
    else:
        atoms_phosp = DNA_PHOSP
        atoms_sugar = DNA_SUGAR
        atoms_nbase = DNA_NBASE
        resnames = ("DT", "DC", "DA", "DG")

    with open(path_json_in, 'r') as file:
        data = json.load(file)

    out = {"phosp": {}, "sugar": {}, "nbase": {}}
    for resname in resnames:
        charges = {k:v[0] for k,v in data[resname].items()}

        out["phosp"][resname] = sum(charges[atomname] for atomname in atoms_phosp)
        out["sugar"][resname] = sum(charges[atomname] for atomname in atoms_sugar)
        out["nbase"][resname] = sum(charges[atomname] for atomname in atoms_nbase[resname])

    with open(path_json_out, 'w') as file:
        json.dump(out, file)


# ------------------------------------------------------------------------------
def main():
    for path_csv in PATHS_CSV:
        csv_to_json(path_csv, path_csv.with_suffix(".json"))

    save_coarse_charges(
        FOLDER_AMBER_CHARGES / "rna.json",
        FOLDER_AMBER_CHARGES / "coarse_rna.json",
        is_rna = True
    )
    save_coarse_charges(
        FOLDER_AMBER_CHARGES / "dna.json",
        FOLDER_AMBER_CHARGES / "coarse_dna.json",
        is_rna = False
    )


################################################################################
if __name__ == "__main__":
    FOLDER_AMBER_CHARGES = Path("volgrids/_data/amber_charges")
    PATHS_CSV = [
        FOLDER_AMBER_CHARGES / name_csv
        for name_csv in ("dna.csv", "rna.csv", "prot.csv")
    ]
    main()


################################################################################
# python3 scripts/dev/get_amber_charges_json.py
