import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class PocketScoreCalculator:
    def __init__(self):
        self.data_pdb = []
        self.data_scores = defaultdict(list)


    # --------------------------------------------------------------------------
    def run(self, folder_data: Path):
        for path_cmap in folder_data.glob("*.cmap"):
            name_pdb = path_cmap.stem
            print(f"Processing {name_pdb}...")

            key_mask = f"{name_pdb}.mask"
            self.data_pdb.append(name_pdb)

            ### boolean grid, points in space that are part of the pocket are "True"
            pocket = np.logical_not(
                vg.GridIO.read_cmap(path_cmap, key_mask).grid.astype(bool)
            )

            keys = set(vg.GridIO.get_cmap_keys(path_cmap)) - {key_mask}
            for key in keys:
                smif = vg.GridIO.read_cmap(path_cmap, key).grid
                kind = key.split('.')[-1]

                if kind == "apbs":
                    neg = smif.copy()
                    pos = smif.copy()
                    neg[neg > 0] = 0
                    pos[pos < 0] = 0
                    self._assign_score("apbs-neg", neg, pocket)
                    self._assign_score("apbs-pos", pos, pocket)

                else:
                    self._assign_score(kind, smif, pocket)


    # --------------------------------------------------------------------------
    def save(self, path_csv: Path):
        print(f"Saving scores to {path_csv}...")
        pd.DataFrame({
            "pdb": self.data_pdb,
            **{kind: scores for kind,scores in self.data_scores.items()}
        }).to_csv(path_csv, index = False)


    # --------------------------------------------------------------------------
    def _assign_score(self, smif_kind: str, smif: vg.Grid, pocket: vg.Grid) -> None:
        volume = len(pocket[pocket])    # the volume is given by the number of points in the pocket
        sum_smif = np.sum(smif[pocket]) # the sum of the smif values in the pocket
        if volume == 0:
            print("Warning: Pocket is empty, no points in the pocket.")
            return


        ##### Option 0: SMIF integral without normalization
        score = np.abs(sum_smif)

        ##### Option 1: SMIF integral is only normalized by the volume of the pocket
        # score = np.abs(sum_smif) / volume

        ##### Option 2: SMIF integral is normalized by its maximum absolute value and by the volume of the pocket
        # maximum_abs = np.max(np.abs(smif))
        # score = np.abs(sum_smif) / (maximum_abs * volume) if maximum_abs > 0 else 0


        self.data_scores[smif_kind].append(score)


################################################################################
if __name__ == "__main__":
    # Run tests/smiffer/pocket_sphere.py before running this script
    # Run this script from the root folder of the repository

    FOLDER_DATA  = Path("testdata/smiffer/pocket_sphere")
    PATH_CSV_OUT = Path("scores.csv")

    psc = PocketScoreCalculator()
    psc.run(FOLDER_DATA)
    psc.save(PATH_CSV_OUT)


################################################################################
