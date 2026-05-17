import sys
import numpy as np
from pathlib import Path
from collections import deque

import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
class ChemTableLigand:
    MAX_BOND_LENGTH = 2.1
    THRESHOLD_PLANARITY = 0.060

    # ------------------------------------------------------------------------------
    def __init__(self, path_pdb_ligand: Path):
        u = mda.Universe(str(path_pdb_ligand))
        self.coords: np.ndarray = u.select_atoms("not name H*").positions
        self.bonds : np.ndarray = self._get_bonds_matrix()
        self.idxs  : np.ndarray = np.arange(len(self.bonds))


    # ------------------------------------------------------------------------------
    def find_flat_rings(self):
        gen_cycles = map(self._find_cycles_dfs, range(len(self.coords)))
        cycles = {
            tuple(sorted(cycle)) : cycle
            for cycle in gen_cycles if cycle is not None
        }

        for cycle in cycles.values():
            is_flat, dev = self._is_flat(cycle)
            msg = "FLAT" if is_flat else "...."
            print(msg, cycle, dev)


    # ------------------------------------------------------------------------------
    def _get_bonds_matrix(self) -> np.ndarray:
        mat0 = self.coords.reshape(1,-1,3)
        mat1 = self.coords.reshape(-1,1,3)
        dist = np.linalg.norm(mat1 - mat0, axis = -1)
        return (0 < dist) & (dist < self.MAX_BOND_LENGTH)


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
    ctl = ChemTableLigand(PATH_PDB_LIGAND)
    ctl.find_flat_rings()


################################################################################
if __name__ == "__main__":
    PATH_PDB_LIGAND = Path(sys.argv[1])
    main()


################################################################################
### check:
# 1raw_ligand_AMP_1
# 6c8e_ligand_EQ1_1
# 6u8u_ligand_Q1V_1
# 8swo_ligand_G3A_2
