import numpy as np
from pathlib import Path
from collections import deque

import volgrids as vg
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class ChemTableLigand:
    ### cycles whose normalized plane normals have a standard deviation
    ### below this value are considered flat enough for stacking
    THRESHOLD_PLANARITY = 0.15

    # ------------------------------------------------------------------------------
    def __init__(self, path_pdb_ligand: Path):
        import MDAnalysis as mda

        u: mda.Universe = vg.Utils.create_mda_universe_quiet(path_pdb_ligand)
        u.guess_TopologyAttrs(to_guess = ["bonds"])
        self.atoms = u.select_atoms("not (name H*) and not water") # [TODO] will presence of ions be an issue?
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
            if is_flat: flat_rings.append(cycle)

            if not su.DEBUG_CHEMTABLE_LIGAND: continue
            msg = "FLAT" if is_flat else "...."
            print(msg, cycle, dev)

        return flat_rings


    # ------------------------------------------------------------------------------
    def gen_chemtable(self) -> str:
        rings = self.find_flat_rings()
        resname = self.atoms[0].resname
        chemtable = '\n'.join((
            "[SELECTION_QUERY]",
            f"resname {resname}", "",
            "[STACKING]",
            ### separate the cycles with extra spaces for readability
            f"{resname}: {'   '.join('-'.join(self.atoms[ring].names) for ring in rings)}",
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

        ### use lists instead of sets to guarantee order of cycle nodes
        ### this is important to have normals consistently oriented when checking planarity
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
    def _is_flat(self, cycle: list[int]) -> tuple[bool, float]:
        normals = []
        for idx,node in enumerate(cycle):
            neighs = self._get_neighs(node)

            ##### immediately reject planarity for tetrahedral atoms
            if len(neighs) > 3:
                return False, np.inf

            ##### add normals for the internal plane of the ring
            i = cycle[(idx-1) % len(cycle)]
            j = node
            k = cycle[(idx+1) % len(cycle)]
            normals.append(self._get_normal(i, j, k))

            if len(neighs) != 3: continue

            ##### add normals for planes defined by any additional substituents
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
