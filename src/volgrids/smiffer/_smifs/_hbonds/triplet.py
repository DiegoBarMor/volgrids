import numpy as np
import MDAnalysis as mda
from MDAnalysis.core.groups import Residue

import volgrids as vg

# ------------------------------------------------------------------------------
def _safe_return_coords(atoms: mda.AtomGroup, sel_string: str):
    sel_atoms = atoms.select_atoms(sel_string)
    if len(sel_atoms) == 0: return None
    return sel_atoms.center_of_geometry()


# //////////////////////////////////////////////////////////////////////////////
class Triplet:
    def __init__(self, res: Residue, interactor: str, tail: str, head: str):
        parts = tail.split('.')
        t0 = parts[0]
        t1 = parts[1] if len(parts) == 2 else ''

        self._t0 = t0
        self._t1 = t1
        self._head = head
        self.interactor = interactor

        self.pos_tail: np.ndarray | None = None
        self.pos_head: np.ndarray | None = None
        self.pos_interactor: np.ndarray | None = None

        self.resname = res.resname.upper()
        self.str_prev_res = f"segid {res.segid} and resid {res.resid - 1}"
        self.str_this_res = f"segid {res.segid} and resid {res.resid    }"
        self.str_next_res = f"segid {res.segid} and resid {res.resid + 1}"

    # --------------------------------------------------------------------------
    def set_pos_tail(self, atoms: mda.AtomGroup) -> np.ndarray | None:
        self.pos_tail = _safe_return_coords(
            atoms, f"name {self._t0} {self._t1}"
        )


    # --------------------------------------------------------------------------
    def set_pos_head(self, atoms: mda.AtomGroup) -> np.ndarray | None:
        self.pos_head = _safe_return_coords(
            atoms, f"name {self._head}"
        )


    # --------------------------------------------------------------------------
    def set_pos_interactor(self, atoms: mda.AtomGroup) -> np.ndarray | None:
        self.pos_interactor = _safe_return_coords(
            atoms, f"name {self.interactor}"
        )


    # --------------------------------------------------------------------------
    def set_pos_tail_custom(self,
        atoms: mda.AtomGroup, query_t0: str, query_t1: str
    ) -> np.ndarray | None:
        self.pos_tail = _safe_return_coords(
            atoms,
            f"({query_t0} and name {self._t0}) or " +\
            f"({query_t1} and name {self._t1})"
        )


    # ------------------------------------------------------------------------------
    def get_interactor_bonded_hydrogens(self, atoms: mda.AtomGroup) -> tuple:
        sel_atoms = atoms.select_atoms(f"name {self.interactor}")
        if len(sel_atoms) == 0:
            return []
        bonded_atoms = [
            (bond.atoms[0] if bond.atoms[0].name != self.interactor else bond.atoms[1])
            for bond in sel_atoms.bonds
        ]
        return tuple(filter(lambda a: a.type == 'H', bonded_atoms))


    # --------------------------------------------------------------------------
    def get_direction_vector(self):
        if (self.pos_tail is None) or (self.pos_head is None):
            return None
        return vg.Math.normalize(self.pos_head - self.pos_tail)


# //////////////////////////////////////////////////////////////////////////////
