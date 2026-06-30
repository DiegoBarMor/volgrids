import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class ParserChemTable:
    def __init__(self, path_table):
        self.resnames: list[str] = []
        self._parser_ini = vg.ParserIni.from_file(path_table)
        self._atoms_hphob: dict[str, dict[str, float]] = {}
        self._names_stk: dict[str, list[str]] = {}
        self._names_hba: dict[str, list[tuple[str, str, str, bool]]] = {}
        self._names_hbd: dict[str, list[tuple[str, str, str, bool]]] = {}
        self._parse_table()


    # --------------------------------------------------------------------------
    def get_atom_hphob(self, atom):
        dict_resid = self._atoms_hphob.get(atom.resname)
        if dict_resid is None: return None
        return dict_resid.get(atom.name)


    # --------------------------------------------------------------------------
    def get_names_stacking(self, resname: str) -> list[str] | None:
        return self._names_stk.get(resname)


    # --------------------------------------------------------------------------
    def get_names_hba(self, resname: str):
        return self._names_hba.get(resname)


    # --------------------------------------------------------------------------
    def get_names_hbd(self, resname: str):
        return self._names_hbd.get(resname)


    # --------------------------------------------------------------------------
    def parse_atom_hphobicity(self, data_ini: vg.ParserIni):
        for resname, str_groups in data_ini.iter_splitted_lines("HYDROPHOBICITY", sep = ':'):
            self._atoms_hphob[resname] = {}
            for group in str_groups.split():
                str_atoms, value = group.split('=')
                value = float(value)
                self._atoms_hphob[resname].update(**{
                    atom.strip() : value for atom in str_atoms.split(',')
                })


    # --------------------------------------------------------------------------
    def parse_names_stacking(self, data_ini: vg.ParserIni):
        for resname, str_cycles in data_ini.iter_splitted_lines("STACKING", sep = ':'):
            self._names_stk[resname] = [
                cycle.replace('-', ' ') for cycle in str_cycles.split()
            ]


    # --------------------------------------------------------------------------
    def parse_names_hbacceptors(self, data_ini: vg.ParserIni):
        for resname, str_triplets in data_ini.iter_splitted_lines("HBACCEPTORS", sep = ':'):
            triplets = map(self._parse_atoms_triplet, str_triplets.split())
            self._names_hba[resname] = [(hba,tail,head,False) for hba,tail,head,_ in triplets] # hbond_fixed must always be False for HBAcceptors


    # --------------------------------------------------------------------------
    def parse_names_hbdonors(self, data_ini: vg.ParserIni):
        for resname, str_triplets in data_ini.iter_splitted_lines("HBDONORS", sep = ':'):
            triplets = list(map(self._parse_atoms_triplet, str_triplets.split()))
            self._names_hbd[resname] = triplets


    # --------------------------------------------------------------------------
    def _parse_table(self):
        """
        Populate the fields of the ParserChemTable instance by parsing the lines of the .chem table file.
        Should only be called once during initialization.
        """

        ### extract values from the lines
        lst_resnames = self._parser_ini.get("RESIDUE_NAMES")
        if lst_resnames is None: raise ValueError("No selection query found in the table file.")
        self.resnames = lst_resnames[0].split()

        self.parse_atom_hphobicity  (self._parser_ini)
        self.parse_names_stacking   (self._parser_ini)
        self.parse_names_hbacceptors(self._parser_ini)
        self.parse_names_hbdonors   (self._parser_ini)


    # --------------------------------------------------------------------------
    @staticmethod
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
