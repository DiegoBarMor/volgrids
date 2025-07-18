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
    def __init__(self, tail: str, head: str, interactor: str):
        parts = tail.split('.')
        t0 = parts[0]
        t1 = parts[1] if len(parts) == 2 else ''

        self._t0 = t0
        self._t1 = t1
        self._head = head
        self._interactor = interactor

        self.pos_tail: np.ndarray = None
        self.pos_head: np.ndarray = None
        self.pos_interactor: np.ndarray = None


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
            atoms, f"name {self._interactor}"
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
    def interactor_is(self, name: str) -> bool:
        return self._interactor == name


    # --------------------------------------------------------------------------
    def get_direction_vector(self):
        if (self.pos_tail is None) or (self.pos_head is None):
            return None
        return vg.Math.normalize(self.pos_head - self.pos_tail)


# //////////////////////////////////////////////////////////////////////////////
