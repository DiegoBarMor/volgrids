import sys
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt

# ------------------------------------------------------------------------------
def main():
    out = "pdb,min,max,avg\n"
    lines = PATH_RAW.read_text().splitlines()
    pairs = list(zip(lines[0::2], lines[1::2]))
    for pair in pairs:
        _,path_grid = pair[0].split(":")
        name_grid = Path(path_grid.strip()).stem
        name_grid = name_grid.replace(".apbs", "")

        raw_stats = pair[1].lstrip(". ").split("; ")
        minval = float(raw_stats[0].split(": ")[1])
        maxval = float(raw_stats[1].split(": ")[1])
        avgval = float(raw_stats[2].split(": ")[1])
        out += f"{name_grid},{minval},{maxval},{avgval}\n"

    PATH_CSV.write_text(out)

    df = pd.read_csv(PATH_CSV)
    df.sort_values("min", inplace = True)
    df.to_csv(PATH_CSV, index = False)

    plt.hist(df["min"], bins = 20)
    plt.show()


################################################################################
if __name__ == "__main__":
    PATH_RAW = Path(sys.argv[1])
    PATH_CSV = Path(sys.argv[2])
    main()


################################################################################
# python3 examples/apbs_summary/parse.py examples/apbs_summary/raw.txt examples/apbs_summary/table.csv
