import sys
import numpy as np
from pathlib import Path
from collections import deque

import MDAnalysis as mda

# //////////////////////////////////////////////////////////////////////////////
class ChemTableLigand:
    MAX_BOND_LENGTH = 2.1

    # ------------------------------------------------------------------------------
    def __init__(self, path_pdb_ligand: Path):
        u = mda.Universe(str(path_pdb_ligand))
        self.coords: np.ndarray = u.select_atoms("not name H*").positions


    # ------------------------------------------------------------------------------
    def find_flat_rings(self):
        mat0 = self.coords.reshape(1,-1,3)
        mat1 = self.coords.reshape(-1,1,3)
        dist = np.linalg.norm(mat1 - mat0, axis = -1)

        bonds = (0 < dist) & (dist < self.MAX_BOND_LENGTH)

        cycles = set(
            self._find_cycles_dfs(bonds, i) for i in range(len(self.coords))
        ) - {None}

        for cycle in cycles:
            if len(cycle) <= 3: continue

            normals = np.zeros((len(cycle), 3))
            for idx,i in enumerate(cycle):
                j = cycle[(idx+1) % len(cycle)]
                k = cycle[(idx+2) % len(cycle)]
                n = np.cross(
                    self.coords[i] - self.coords[j],
                    self.coords[k] - self.coords[j]
                )
                n /= np.linalg.norm(n)
                normals[idx] = n

            dev = np.linalg.norm(
                np.std(normals, axis = 0)
            )
            print(cycle, "DEVIATION:", dev)


    # ------------------------------------------------------------------------------
    @classmethod
    def _find_cycles_dfs(cls, graph: np.ndarray, idx_start: int) -> tuple[int]:
        paths  : deque[set] = deque()
        queue  : deque[int] = deque()
        parents: deque[int] = deque()

        paths.append(set())
        queue.append(idx_start)
        parents.append(-1)

        idxs = np.arange(len(graph))
        while queue:
            node = int(queue.popleft())
            path = paths.popleft()
            parent = parents.popleft()

            neighs = [int(i) for i in idxs[graph[node]] if i != parent]

            if len(path) and node == idx_start:
                return tuple(sorted(path))

            if node in path: continue

            path.add(node)

            queue.extend(neighs)
            paths.extend([path.copy() for _ in neighs])
            parents.extend([node for _ in neighs])


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
