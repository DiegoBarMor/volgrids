import re
from collections import defaultdict

############################### PARSING UTILITIES ##############################
# ------------------------------------------------------------------------------
def _extract_section_between_markers(superstr: str, start_regex: str, end_regex: str) -> str:
    match_start = re.search(start_regex, superstr)
    if match_start is None:
        raise ValueError(f"Start marker '{start_regex}' not found in the string.")

    _,i0 = match_start.span()
    sliced_str = superstr[i0:]

    match_end = re.search(end_regex, sliced_str)
    if match_end is None:
        return sliced_str.strip()

    i1,_ = match_end.span()
    return sliced_str[:i1].strip()


# ------------------------------------------------------------------------------
def _split_line(line: str) -> tuple[str, str]:
    line = line.split('#')[0].strip() # Remove comments
    pair = tuple(map(str.strip, line.split(':')))
    if len(pair) != 2:
        raise ValueError(f"Line '{line}' does not contain ':'")
    return pair


# ------------------------------------------------------------------------------
def _is_empty_line(line: str) -> bool:
    return not line.strip() or line.startswith('#')


# ------------------------------------------------------------------------------
def _parse_atoms_triplet(triplet: str) -> tuple[str, str, str]:
    try:
        str_antecedents, interactor = triplet.split('-')
        antecedents = str_antecedents.split('.')
        a0 = antecedents[0].strip()
        a1 = antecedents[1].strip() if len(antecedents) > 1 else ''
        return a0, a1, interactor.strip()
    except ValueError as e:
        raise ValueError(f"Triplet '{triplet}' is not in the expected formats 'A-X' or 'A.B-X'") from e


# //////////////////////////////////////////////////////////////////////////////
class ChemTable:
    def __init__(self, path_table):
        self.selection_query: str = ''
        self._residues_hphob: dict[str, float] = {}
        self._atoms_hphob: defaultdict[str, dict[str, float]] = defaultdict(dict)
        self._names_stk: dict[str, str] = {}
        self._names_hba: dict[str, list[tuple[str, str, str]]] = {}
        self._names_hbd: dict[str, list[tuple[str, str, str]]] = {}
        self._parse_table(path_table)


    # --------------------------------------------------------------------------
    def get_residue_hphob(self, atom):
        return self._residues_hphob.get(atom.resname)


    # --------------------------------------------------------------------------
    def get_atom_hphob(self, atom):
        dict_resid = self._atoms_hphob.get(atom.resname)
        if dict_resid is None: return None
        return dict_resid.get(atom.name)


    # --------------------------------------------------------------------------
    def get_names_stacking(self, resname: str):
        return self._names_stk.get(resname)


    # --------------------------------------------------------------------------
    def get_names_hba(self, resname: str):
        return self._names_hba.get(resname)


    # --------------------------------------------------------------------------
    def get_names_hbd(self, resname: str):
        return self._names_hbd.get(resname)


    # --------------------------------------------------------------------------
    def _parse_table(self, path_table):
        with open(path_table, 'r') as file:
            raw_table = file.read()

        ### separate raw string into sections
        raw_selection_query = _extract_section_between_markers(
            raw_table, r"\[SELECTION_QUERY\]", r"\[\w*\]"
        )
        raw_residues_hphob = _extract_section_between_markers(
            raw_table, r"\[RES_HPHOBICITY\]", r"\[\w*\]"
        )
        raw_atoms_hphob = _extract_section_between_markers(
            raw_table, r"\[ATOM_HPHOBICITY\]", r"\[\w*\]"
        )
        raw_names_stk = _extract_section_between_markers(
            raw_table, r"\[NAMES_STACKING\]", r"\[\w*\]"
        )
        raw_names_hba = _extract_section_between_markers(
            raw_table, r"\[NAMES_HBACCEPTORS\]", r"\[\w*\]"
        )
        raw_names_hbd = _extract_section_between_markers(
            raw_table, r"\[NAMES_HBDONORS\]", r"\[\w*\]"
        )

        ### extract values from the lines
        self.selection_query = raw_selection_query.strip()

        for line in raw_residues_hphob.splitlines():
            if _is_empty_line(line): continue
            resname, value = _split_line(line)
            self._residues_hphob[resname] = float(value)

        for line in raw_atoms_hphob.splitlines():
            if _is_empty_line(line): continue
            names, value = _split_line(line)
            resname, atomname = names.split('/')
            self._atoms_hphob[resname][atomname] = float(value)

        for line in raw_names_stk.splitlines():
            if _is_empty_line(line): continue
            resname, atomnames = _split_line(line)
            self._names_stk[resname] = atomnames

        for line in raw_names_hba.splitlines():
            if _is_empty_line(line): continue
            resname, str_triplets = _split_line(line)
            triplets = list(map(_parse_atoms_triplet, str_triplets.split()))
            self._names_hba[resname] = triplets

        for line in raw_names_hbd.splitlines():
            if _is_empty_line(line): continue
            resname, str_triplets = _split_line(line)
            triplets = list(map(_parse_atoms_triplet, str_triplets.split()))
            self._names_hbd[resname] = triplets

        ### expand the atoms declared with the '*' wildcard to all residues (hydrophobicity)
        hphob_wildcard = self._atoms_hphob.get('*')
        if hphob_wildcard is None: return

        del self._atoms_hphob['*']
        for res_dict in self._atoms_hphob.values():
            res_dict.update(hphob_wildcard)


# //////////////////////////////////////////////////////////////////////////////
