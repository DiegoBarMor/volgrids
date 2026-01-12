from collections import defaultdict

import volgrids as vg

# ------------------------------------------------------------------------------
def _parse_atoms_triplet(triplet: str) -> tuple[str, str, str, str, bool]:
    def _assert(condition: bool):
        assert condition, \
            f"Triplet '{triplet}' is not in the expected formats 'I=T->H' or 'I=T0.T1->H'."

    stripped = triplet.strip('!')
    hbond_fixed = stripped != triplet

    parts = stripped.split('=')
    _assert(len(parts) == 2)

    interactor, direction = parts
    parts = direction.split("->")
    _assert(len(parts) == 2)

    tail, head = parts
    _assert(tail and head and interactor)

    ### valid syntax?  | yes | no |
    ### 'i=t->h'       |  X  |    |
    ### 'i=t0.t1->h'   |  X  |    |
    ### 'i=t0.t1.t2->h'|  X  |    |
    ### 'i=t0.->h'     |     | X  |
    ### 'i=.t1->h'     |     | X  |
    tail_points = tail.split('.')
    _assert(len(tail_points) > 0 and all(tail_points))

    return interactor, tail_points, head, hbond_fixed


# //////////////////////////////////////////////////////////////////////////////
class ParserChemTable:
    def __init__(self, path_table):
        self.selection_query: str = ''
        self._parser_ini = vg.ParserIni(path_table)
        self._residues_hphob: dict[str, float] = {}
        self._atoms_hphob: defaultdict[str, dict[str, float]] = defaultdict(dict)
        self._names_stk: dict[str, str] = {}
        self._names_hba: dict[str, list[tuple[str, str, str, bool]]] = {}
        self._names_hbd: dict[str, list[tuple[str, str, str, bool]]] = {}
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
        query = self._parser_ini.get("SELECTION_QUERY")
        if query is None: raise ValueError("No selection query found in the table file.")
        self.selection_query = query[0]

        for resname, value in self._parser_ini.iter_splitted_lines("RES_HPHOBICITY", sep = ':'):
            self._residues_hphob[resname] = float(value)

        for names, value in self._parser_ini.iter_splitted_lines("ATOM_HPHOBICITY", sep = ':'):
            resname, atomname = names.split('/')
            self._atoms_hphob[resname][atomname] = float(value)

        for resname, atomnames in self._parser_ini.iter_splitted_lines("NAMES_STACKING", sep = ':'):
            self._names_stk[resname] = atomnames

        for resname, str_triplets in self._parser_ini.iter_splitted_lines("NAMES_HBACCEPTORS", sep = ':'):
            triplets = map(_parse_atoms_triplet, str_triplets.split())
            self._names_hba[resname] = [(hba,tail,head,False) for hba,tail,head,_ in triplets] # hbond_fixed must always be False for HBAcceptors

        for resname, str_triplets in self._parser_ini.iter_splitted_lines("NAMES_HBDONORS", sep = ':'):
            triplets = list(map(_parse_atoms_triplet, str_triplets.split()))
            self._names_hbd[resname] = triplets

        ### expand the atoms declared with the '*' wildcard to all residues (hydrophobicity)
        hphob_wildcard = self._atoms_hphob.get('*')
        if hphob_wildcard is None: return

        del self._atoms_hphob['*']
        for res_dict in self._atoms_hphob.values():
            res_dict.update(hphob_wildcard)


# //////////////////////////////////////////////////////////////////////////////
