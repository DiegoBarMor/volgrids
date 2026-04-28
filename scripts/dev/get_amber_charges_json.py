import json
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------------------------
def main():
    for path_csv in PATHS_CSV:
        df = pd.read_csv(path_csv)
        resnames = df["resname"].unique()
        out = {}
        for resname in resnames:
            df_res = df[df["resname"] == resname]
            out[resname] = dict(zip(
                df_res["atomname"],
                list(zip(df_res["charge"], df_res["radius"]))
            ))

        with open(path_csv.with_suffix(".json"), 'w') as file:
            json.dump(out, file, indent = 4)


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
