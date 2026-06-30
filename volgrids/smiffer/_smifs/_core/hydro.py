from abc import ABC

from .smif import Smif

# //////////////////////////////////////////////////////////////////////////////
class SmifHydro(Smif, ABC):
    def iter_particles(self):
        for atom in self.mm.get_atoms_insphere():
            factor_atom = self.mm.chemtable.get_atom_hphob(atom)
            if factor_atom is None: continue # skip unknown atoms
            yield atom, factor_atom #/ len(atom.residue.atoms)


# //////////////////////////////////////////////////////////////////////////////
