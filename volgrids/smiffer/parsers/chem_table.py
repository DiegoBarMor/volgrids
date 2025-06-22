import volgrids as vg
from collections import defaultdict

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
class ChemTable():
    def __init__(self, path_table):
        self.ini_parser = vg.IniParser(path_table)
        self.selection_query: str = ''
        self._residues_hphob: dict[str, float] = {}
        self._atoms_hphob: defaultdict[str, dict[str, float]] = defaultdict(dict)
        self._names_stk: dict[str, str] = {}
        self._names_hba: dict[str, list[tuple[str, str, str]]] = {}
        self._names_hbd: dict[str, list[tuple[str, str, str]]] = {}
        self._parse_table()


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
    def _parse_table(self):
        ### extract values from the lines
        query = self.ini_parser.get("SELECTION_QUERY")
        if query is None: raise ValueError("No selection query found in the table file.")
        self.selection_query = query[0]

        for resname, value in self.ini_parser.iter_splitted_lines("RES_HPHOBICITY", sep = ':'):
            self._residues_hphob[resname] = float(value)

        for names, value in self.ini_parser.iter_splitted_lines("ATOM_HPHOBICITY", sep = ':'):
            resname, atomname = names.split('/')
            self._atoms_hphob[resname][atomname] = float(value)

        for resname, atomnames in self.ini_parser.iter_splitted_lines("NAMES_STACKING", sep = ':'):
            self._names_stk[resname] = atomnames

        for resname, str_triplets in self.ini_parser.iter_splitted_lines("NAMES_HBACCEPTORS", sep = ':'):
            triplets = list(map(_parse_atoms_triplet, str_triplets.split()))
            self._names_hba[resname] = triplets

        for resname, str_triplets in self.ini_parser.iter_splitted_lines("NAMES_HBDONORS", sep = ':'):
            triplets = list(map(_parse_atoms_triplet, str_triplets.split()))
            self._names_hbd[resname] = triplets

        ### expand the atoms declared with the '*' wildcard to all residues (hydrophobicity)
        hphob_wildcard = self._atoms_hphob.get('*')
        if hphob_wildcard is None: return

        del self._atoms_hphob['*']
        for res_dict in self._atoms_hphob.values():
            res_dict.update(hphob_wildcard)


# //////////////////////////////////////////////////////////////////////////////
