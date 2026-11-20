import numpy as np

import volgrids as vg
import volgrids._vendors.molsimple as ms

# //////////////////////////////////////////////////////////////////////////////
class Triplet:
    def __init__(self,
        interactor: str,
        tail: tuple[str],
        head: str,
        hbond_fixed: bool,
    ):
        self._tail = tail
        self._head = head
        self.interactor = interactor
        self.hbond_fixed = hbond_fixed

        self.pos_tail: np.ndarray | None = None
        self.pos_head: np.ndarray | None = None
        self.pos_interactor: np.ndarray | None = None

        self.resname: str = ""
        self.residue_prev: ms.ParticleGroup | None = None
        self.residue_this: ms.ParticleGroup | None = None
        self.residue_next: ms.ParticleGroup | None = None


    # --------------------------------------------------------------------------
    def has_prev_res(self) -> bool:
        return self.residue_prev is not None


    # --------------------------------------------------------------------------
    def has_next_res(self) -> bool:
        return self.residue_next is not None


    # --------------------------------------------------------------------------
    def reset_posits(self) -> None:
        self.pos_tail = None
        self.pos_head = None
        self.pos_interactor = None


    # --------------------------------------------------------------------------
    def set_pos_tail(self, tail_override: tuple[str] = None) -> None:
        tail = self._tail if tail_override is None else tail_override
        self.pos_tail = self._safe_cog(self.residue_this, *tail)


    # --------------------------------------------------------------------------
    def set_pos_head(self) -> None:
        self.pos_head = self._safe_cog(self.residue_this, self._head)


    # --------------------------------------------------------------------------
    def set_pos_interactor(self) -> None:
        self.pos_interactor = self._safe_cog(self.residue_this, self.interactor)


    # --------------------------------------------------------------------------
    def set_pos_tail_custom(self,
        res_0: ms.ParticleGroup, res_1: ms.ParticleGroup
    ) -> np.ndarray | None:
        """
        Reserved for special cases in protein and RNA.
        """
        assert len(self._tail) == 2, "Custom tail position can only be set for triplets with two tail points."
        t0, t1 = self._tail
        cog_0 = self._safe_cog(res_0, t0)
        cog_1 = self._safe_cog(res_1, t1)
        self.pos_tail = (cog_0 + cog_1) / 2 if (cog_0 is not None and cog_1 is not None) else None


    # --------------------------------------------------------------------------
    def get_direction_vector(self):
        if (self.pos_tail is None) or (self.pos_head is None):
            return None
        return vg.Math.normalize(self.pos_head - self.pos_tail)


    # ------------------------------------------------------------------------------
    @staticmethod
    def _safe_cog(residue: ms.ParticleGroup|None, *atom_names: str) -> np.ndarray | None:
        if residue is None: return
        atoms = residue.select_name(*atom_names)
        if len(atoms): return atoms.get_center_of_geometry()


# //////////////////////////////////////////////////////////////////////////////
