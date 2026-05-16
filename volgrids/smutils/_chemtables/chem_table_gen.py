from rdkit import Chem
from rdkit.Chem import Draw

# //////////////////////////////////////////////////////////////////////////////
class ChemTableGenerator:
    def __init__(self, molecule: Chem.Mol):
        self.molecule = molecule

    # --------------------------------------------------------------------------
    @classmethod
    def from_smiles(cls, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        return cls(mol)

    # --------------------------------------------------------------------------
    @classmethod
    def from_pdb(cls, path_pdb: str):
        mol = Chem.MolFromPDBFile(path_pdb, removeHs = False)
        return cls(mol)

    # --------------------------------------------------------------------------
    def show_2d(self):
        img = Draw.MolToImage(self.molecule)
        img.show()


# //////////////////////////////////////////////////////////////////////////////



OCTET = 8
ELECTRONS_PER_BOND = 2
ELECTRONS_PER_LONE_PAIR = 2
HBOND_HETEROATOMS = set(('O', 'N', 'P', 'S'))

def infer_lone_pairs(atom: Chem.Atom) -> int:
    # nbonds_explicit = atom.GetValence(Chem.ValenceType.EXPLICIT) # bonds to other atoms
    # nbonds_implicit = atom.GetValence(Chem.ValenceType.IMPLICIT) # hydrogens
    # ncovalent_electrons = (nbonds_explicit + nbonds_implicit) * ELECTRONS_PER_BOND
    nbonds = atom.GetTotalValence()
    nlone_pairs = (OCTET - nbonds*ELECTRONS_PER_BOND) // ELECTRONS_PER_LONE_PAIR
    return nlone_pairs


###### SOME WIP
# def get_nbonded_hydrogens(atom: Chem.Atom) -> int:
#     nh = 0
#     for a1 in atom.GetNeighbors():
#         if a1.GetSymbol() == 'H':
#             nh += 1
#     return nh

# path_cif = "testdata/_input/ligands/mid_size/BEB.cif"
# mol = Chem.MolFromPDBFile("testdata/_input/ligands/mid_size/beb.pdb", removeHs = False)
# mol = Chem.MolFromPDBFile("testdata/smiffer/ligands/all_interactions/tsc.pdb", removeHs = False)



############ SOME MORE TESTING
# smiles_thiamine = "OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N"
# smiles_beb = "O=C(NC2c1ccccc1CC2O)C(OCc3ccccc3)C(O)C(O)C(OCc4ccccc4)C(=O)NC6c5ccccc5CC6O"
# smiles_salt = "[Cu+2].[O-]S(=O)(=O)[O-]"

# generator = ChemTableGenerator.from_smiles(smiles_salt)

# # for atom in mol.GetAtoms():
# #     atom: Chem.Atom
# #     if atom.GetSymbol().upper() not in HBOND_HETEROATOMS: continue

# #     nlone_pairs = infer_lone_pairs(atom)
# #     if nlone_pairs <= 0: continue


# #     print(atom.GetSymbol(), nlone_pairs)

# #     # for a1 in atom.GetNeighbors():
# #     #     print(a1.GetIdx(), a1.GetSymbol(), end = ' ')
# #     # exit()


# # # for bond in mol.GetBonds():
# # #     bond: Chem.Bond
# # #     print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType(), bond.GetIsAromatic())



######################## UI WITH THE OLD SCHEMA
# # --------------------------------------------------------------------------
# def assign_globals(self):
#     self._set_help_str(
#         "usage: python3 run/utils.py [mode] [options...]",
#         "Available modes:",
#         "  chemgen - Generate a SMIFer chem table for an arbitrary molecule (e.g. ligands).",
#         "Run 'python3 run/utils.py [mode] --help' for more details on each mode.",
#         "Running 'python3 run/utils.py' without a valid mode will display this help message.",
#     )
#     if self._has_param_kwds("help") and not self._has_params_pos():
#         self._exit_with_help(0)

#     vu.MODE = self._safe_get_param_pos(0).lower()
#     func: callable = self._safe_map_value(vu.MODE,
#         chemgen = self._parse_chemgen,
#     )
#     func()


# # --------------------------------------------------------------------------
# def _parse_chemgen(self) -> None:
#     self._set_help_str(
#         "usage: python3 run/utils.py chemgen [options...]",
#         "Available options:",
#         "-h, --help    Show this help message and exit.",
#         "-p, --pdb     Path to a PDB file with the isolated molecule structure. Use either this or the '--smiles' option.",
#         "-s, --smiles  SMILES string of the molecule. Use either this or the '--pdb' option.",
#         "-o, --output  Path to the folder where the output CHEM tables should be stored. If not provided, the parent folder of the input PDB structure file will be used (mandatory if using SMILES input).",
#     )
#     if self._has_param_kwds("help"):
#         self._exit_with_help(0)

#     flag_pdb    = self._has_param_kwds("pdb")
#     flag_smiles = self._has_param_kwds("smiles")

#     if flag_pdb == flag_smiles:
#         self._exit_with_help(-1, "You must provide either the '--pdb' or the '--smiles' option (but not both).")

#     if flag_pdb:
#         vu.INPUT_KIND = "pdb"
#         vu.PATH_PDB = self._safe_path_file_in(
#             self._safe_get_param_kwd("pdb"),
#             err_msg = "No valid PDB file provided. Provide a path to the PDB file using the '--pdb' option."
#         )

#     if flag_smiles:
#         vu.INPUT_KIND = "smiles"
#         vu.STR_SMILES = self._safe_get_param_kwd("smiles")

#     if self._has_param_kwds("output"):
#         vu.PATH_OUTPUT = self._safe_path_folder_out(
#             self._safe_get_param_kwd("output")
#         )
#     elif flag_pdb:
#         vu.PATH_OUTPUT = vu.PATH_PDB.parent
#     else:
#         self._exit_with_help(-1, "You must provide the '--output' option when using SMILES input.")
