import time

from src.utilities import save_metadata
from src.args_manager import process_args
from src.mol_systems import MS_Whole, MS_PocketSphere
from src.grids.potential_grids import PotentialGrid
from src.grids.apbs import PPG_APBS
from src.grids.hbonds import SPG_HB_Accepts, SPG_HB_Donors
from src.grids.hydro import SPG_Hydrophobic, SPG_Hydrophilic
from src.grids.stacking import SPG_Stacking

from settings import TRIMMING_DIST_LARGE, TRIMMING_DIST_SMALL


################################################################################
if __name__ == "__main__":
    args = process_args()
    metadata = vars(args)
    metadata["meta"] = args.out / f"{args.pdb.stem}.meta.json"
    save_metadata(metadata)

    print(f">>> Now processing '{args.pdb.stem}' ({'nucleic' if args.rna else 'protein'})...", flush = True)

    ms_class = MS_Whole if metadata["whole"] else MS_PocketSphere

    # start = time.time()
    ms_trimming_large = ms_class(metadata, TRIMMING_DIST_LARGE)
    ms_trimming_small = ms_class(metadata, TRIMMING_DIST_SMALL)
    pg_stacking  = SPG_Stacking   (ms_trimming_large)
    pg_hba       = SPG_HB_Accepts (ms_trimming_large)
    pg_hbd       = SPG_HB_Donors  (ms_trimming_large)
    pg_hbphob    = SPG_Hydrophobic(ms_trimming_large)
    pg_hphil     = SPG_Hydrophilic(ms_trimming_small)
    pg_apbs      = PPG_APBS       (ms_trimming_large)
    pg_hydrodiff = PotentialGrid.grid_diff(pg_hbphob, pg_hphil, "hydrodiff")

    # print("Elapsed:", time.time() - start, flush = True)


################################################################################
