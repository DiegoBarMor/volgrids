from pathlib import Path

from volgrids._vendors import molutils as mu

# //////////////////////////////////////////////////////////////////////////////
class ResiduesNucleic:
    STACKING_THRESHOLD = 2 # [TODO] config?

    # --------------------------------------------------------------------------
    @classmethod
    def get_residues_bp_canonical(cls, path_pdb: Path) -> set[str]:
        from rnapolis.common import Structure2D

        structure2d: Structure2D = cls._get_rnapolis_struct2d(path_pdb)
        base_pairs = [bp for bp in structure2d.base_pairs if bp.saenger is not None]
        bp_full_names = [
            [bp.nt1.full_name, bp.nt2.full_name] for bp in base_pairs
            if bp.lw.value == "cWW" and bp.saenger.value in ("XIX", "XX")
        ]
        bps = {bp[0] for bp in bp_full_names} | {bp[1] for bp in bp_full_names}
        return {cls._rnapolis_residue(bp) for bp in bps}


    # --------------------------------------------------------------------------
    @classmethod
    def get_residues_stacking(cls, path_pdb: Path) -> set[str]:
        """Return residues that are participating in at least `STACKING_THRESHOLD` stacking interactions."""
        from collections import Counter
        from rnapolis.common import Structure2D

        structure2d: Structure2D = cls._get_rnapolis_struct2d(path_pdb)
        stackings = [bp for bp in structure2d.stackings if bp.topology is not None]
        stk_full_names = [
            [stk.nt1.full_name, stk.nt2.full_name] for stk in stackings
        ]
        resids = [stk[0] for stk in stk_full_names] + [stk[1] for stk in stk_full_names]
        return set(
            cls._rnapolis_residue(residue) for residue, count in Counter(resids).items()
            if count >= cls.STACKING_THRESHOLD
        )


    # --------------------------------------------------------------------------
    @classmethod
    def get_residues_nobp(cls, path_pdb: str | Path) -> str:
        path_pdb = Path(path_pdb)
        resids_all = set(mu.List.residues(path_pdb))
        resids_canonical = cls.get_residues_bp_canonical(path_pdb)
        return ' '.join(sorted(resids_all - resids_canonical))


    # --------------------------------------------------------------------------
    @classmethod
    def get_residues_nostk(cls, path_pdb: str | Path) -> str:
        path_pdb = Path(path_pdb)
        resids_all = set(mu.List.residues(path_pdb))
        resids_stacking = cls.get_residues_stacking(path_pdb)
        return ' '.join(sorted(resids_all - resids_stacking))


    # --------------------------------------------------------------------------
    @staticmethod
    def _rnapolis_residue(s):
        chain,residue = s.split('.')
        while residue:
            if residue[0].isdigit() or residue[0] == '-': break
            residue = residue[1:]
        return f"{chain}.{residue}"


    # --------------------------------------------------------------------------
    @staticmethod
    def _get_rnapolis_struct2d(path_pdb: Path):
        from rnapolis.annotator import  (
            extract_base_interactions,
            handle_input_file,
            read_3d_structure,
        )

        file = handle_input_file(path_pdb)
        structure3d = read_3d_structure(file, None)
        base_interactions = extract_base_interactions(structure3d)
        structure2d, _ = structure3d.extract_secondary_structure(base_interactions)
        return structure2d


# //////////////////////////////////////////////////////////////////////////////
