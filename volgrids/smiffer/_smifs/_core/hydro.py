from abc import ABC

from .smif import Smif

# //////////////////////////////////////////////////////////////////////////////
class SmifHydro(Smif, ABC):
    def iter_particles(self):
        for atom in self.ms.get_relevant_atoms():
            factor_atom = self.ms.chemtable.get_atom_hphob(atom)

            if factor_atom is None: continue # skip atoms with unknown name

            yield atom, factor_atom #/ len(atom.residue.atoms)


# //////////////////////////////////////////////////////////////////////////////
