import io
import logging
import warnings
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import MDAnalysis as mda

from rnapolis.common import Structure2D
from rnapolis.annotator import  (
    extract_base_interactions,
    handle_input_file,
    read_3d_structure,
)

# //////////////////////////////////////////////////////////////////////////////
class RNAResids:
    # --------------------------------------------------------------------------
    @classmethod
    def get_resids_bp_canonical(cls, path_pdb: Path) -> set[int]:
        structure2d = cls._get_rnapolis_struct2d(path_pdb)
        base_pairs = [bp for bp in structure2d.base_pairs if bp.saenger is not None]
        bp_full_names = [
            [bp.nt1.full_name, bp.nt2.full_name] for bp in base_pairs
            if bp.lw.value == "cWW" and bp.saenger.value in ("XIX", "XX")
        ]
        resids_0 = { int(bp[0][3:]) for bp in bp_full_names }
        resids_1 = { int(bp[1][3:]) for bp in bp_full_names }
        return resids_0 | resids_1


    # --------------------------------------------------------------------------
    @classmethod
    def get_all_resids(cls, path_pdb) -> set[int]:
        u = cls._create_mda_universe_quiet(path_pdb)
        return set(int(i) for i in u.residues.resids)


    # --------------------------------------------------------------------------
    @classmethod
    def get_resids_nonbp(cls, path_pdb: str | Path) -> str:
        path_pdb = Path(path_pdb)
        idxs_canonical = cls.get_resids_bp_canonical(path_pdb)
        idxs_all = cls.get_all_resids(path_pdb)
        idxs_nonpb = sorted(idxs_all - idxs_canonical)
        return ' '.join(str(i) for i in idxs_nonpb)


    # --------------------------------------------------------------------------
    @staticmethod
    def _create_mda_universe_quiet(path_pdb) -> mda.Universe:
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


    # --------------------------------------------------------------------------
    @staticmethod
    def _get_rnapolis_struct2d(path_pdb: Path) -> Structure2D:
        file = handle_input_file(path_pdb)
        structure3d = read_3d_structure(file, None)
        base_interactions = extract_base_interactions(structure3d)
        structure2d, _ = structure3d.extract_secondary_structure(
            base_interactions, False, False
        )
        return structure2d


# //////////////////////////////////////////////////////////////////////////////
