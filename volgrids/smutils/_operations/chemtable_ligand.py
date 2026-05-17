import sys
import numpy as np
from pathlib import Path
from collections import deque

import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
class ChemTableLigand:
    THRESHOLD_PLANARITY = 0.15

    # ------------------------------------------------------------------------------
    def __init__(self, path_pdb_ligand: Path):
        u = mda.Universe(str(path_pdb_ligand))
        u.guess_TopologyAttrs(to_guess = ["bonds"])
        self.atoms = u.select_atoms("not (name H*)")
        self.coords: np.ndarray = self.atoms.positions
        self.bonds : np.ndarray = self._get_bonds_matrix()
        self.idxs  : np.ndarray = np.arange(len(self.bonds))


    # ------------------------------------------------------------------------------
    def find_flat_rings(self):
        gen_cycles = map(self._find_cycles_dfs, range(len(self.coords)))
        cycles = {
            tuple(sorted(cycle)) : cycle
            for cycle in gen_cycles if cycle is not None
        }

        flat_rings = []
        for cycle in cycles.values():
            is_flat, dev = self._is_flat(cycle)
            msg = "FLAT" if is_flat else "...."
            print(msg, cycle, dev)
            if is_flat: flat_rings.append(cycle)

        return flat_rings


    # ------------------------------------------------------------------------------
    def gen_chemtable(self) -> str:
        rings = self.find_flat_rings()
        resname = self.atoms[0].resname
        chemtable = '\n'.join((
            "[SELECTION_QUERY]",
            f"resname {resname}", "",
            "[NAMES_STACKING]",
            *(
                f"{resname}: {' '.join(self.atoms[ring].names)}"
                for ring in rings
            )
        ))
        return chemtable


    # ------------------------------------------------------------------------------
    def _get_bonds_matrix(self) -> np.ndarray:
        mat = np.zeros((len(self.coords), len(self.coords)), dtype = bool)
        atomid_to_idx = {int(a.id) : i for i,a in enumerate(self.atoms)}
        for bond in self.atoms.bonds:
            i = atomid_to_idx.get(int(bond.atoms[0].id))
            j = atomid_to_idx.get(int(bond.atoms[1].id))
            if i is None: continue # skip bonds with excluded atoms i.e. hydrogens
            if j is None: continue
            mat[i,j] = True
            mat[j,i] = True
        return mat


    # ------------------------------------------------------------------------------
    def _get_neighs(self, idx: int) -> np.ndarray:
        return self.idxs[self.bonds[idx]]


    # ------------------------------------------------------------------------------
    def _get_normal(self, i: int, j: int, k: int) -> np.ndarray:
        n = np.cross(
            self.coords[i] - self.coords[j],
            self.coords[k] - self.coords[j]
        )
        return n / np.linalg.norm(n)


    # ------------------------------------------------------------------------------
    def _find_cycles_dfs(self, idx_start: int) -> tuple[int,...]:
        paths  : deque[list] = deque()
        queue  : deque[int]  = deque()
        parents: deque[int]  = deque()

        paths.append([])
        queue.append(idx_start)
        parents.append(-1)

        while queue:
            node = int(queue.popleft())
            path = paths.popleft()
            parent = parents.popleft()

            neighs = [int(i) for i in self._get_neighs(node) if i != parent]

            if len(path) and node == idx_start:
                return path

            if node in path: continue

            path.append(node)

            queue.extend(neighs)
            paths.extend([path.copy() for _ in neighs])
            parents.extend([node for _ in neighs])


    # ------------------------------------------------------------------------------
    def _is_flat(self, cycle: list[int]) -> bool:
        normals = []
        for idx,node in enumerate(cycle):
            neighs = self._get_neighs(node)

            if len(neighs) > 3:
                return False, np.inf

            i = cycle[(idx-1) % len(cycle)]
            j = node
            k = cycle[(idx+1) % len(cycle)]
            normals.append(self._get_normal(i, j, k))

            if len(neighs) != 3: continue

            l = neighs[0] if neighs[0] not in (i,k) else (
                neighs[1] if neighs[1] not in (i,k) else
                neighs[2]
            )
            normals.append(self._get_normal(l, j, i))
            normals.append(self._get_normal(k, j, l))

        dev = np.linalg.norm(
            np.std(normals, axis = 0)
        )
        return dev < self.THRESHOLD_PLANARITY, dev


# //////////////////////////////////////////////////////////////////////////////
# ------------------------------------------------------------------------------
def main():
    PATH_CHEM_OUT.write_text(
        ChemTableLigand(PATH_PDB_LIGAND).gen_chemtable()
    )


################################################################################
if __name__ == "__main__":
    PATH_PDB_LIGAND = Path(sys.argv[1])
    PATH_CHEM_OUT = PATH_PDB_LIGAND.with_suffix(".chem")
    main()


################################################################################
### check:
# 1raw_ligand_AMP_1
# 6c8e_ligand_EQ1_1
# 6u8u_ligand_Q1V_1
# 8swo_ligand_G3A_2
