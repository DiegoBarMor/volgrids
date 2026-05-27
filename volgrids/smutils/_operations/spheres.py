import numpy as np
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class Spheres:
    # --------------------------------------------------------------------------
    @staticmethod
    def find_spheres(path_pdb: Path, path_traj: Path, query: str, radius_extra: float) -> str:
        def _sphere_info():
            coords = atoms.positions
            cog = np.mean(coords, axis = 0)
            max_dist = np.max(np.linalg.norm(coords - cog, axis = 1))
            radius = max_dist + radius_extra
            return f"{cog[0]:.3f} {cog[1]:.3f} {cog[2]:.3f} {radius:.3f}"

        u = vg.create_mda_universe_quiet(path_pdb, path_traj)
        atoms = u.select_atoms(query)
        if len(atoms) == 0:
            raise ValueError(f"No atoms found for query: {query}")

        return ' '.join(_sphere_info() for _ in u.trajectory)


# //////////////////////////////////////////////////////////////////////////////
