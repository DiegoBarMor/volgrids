from abc import ABC

import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifHydro(sm.Smif, ABC):
    def iter_particles(self):
        for atom in self.ms.get_relevant_atoms():
            factor_res  = self.ms.chemtable.get_residue_hphob(atom)
            factor_atom = self.ms.chemtable.get_atom_hphob(atom)

            if (factor_res is None) and (factor_atom is None):
                continue # skip atoms with unknown name and resname

            if factor_res  is None: factor_res  = 1
            if factor_atom is None: factor_atom = 1

            yield atom, factor_res * factor_atom #/ len(atom.residue.atoms)


# //////////////////////////////////////////////////////////////////////////////
