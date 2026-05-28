import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smutils as su
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSpheres(vg.AppSubcommand):
    # -------------------------------------------------------------------------- UI SECTION
    def run(self):
        self.app_run_find()


    # --------------------------------------------------------------------------
    def app_run_find(self):
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_traj = self.main.get_arg_path("path_traj",
            assertion = fy.PathAssertion.FILE_IN, allow_none = True
        )
        query = self.main.get_arg_str("query")
        radius_extra = self.main.get_arg_float("radius_extra")

        print(su.AppSpheres.find(path_in, path_traj, query, radius_extra))


    # -------------------------------------------------------------------------- LOGIC SECTION
    @staticmethod
    def find(path_pdb: Path, path_traj: Path|None, query: str, radius_extra: float) -> str:
        def get_sphere_info():
            coords = atoms.positions
            cog = np.mean(coords, axis = 0)
            max_dist = np.max(np.linalg.norm(coords - cog, axis = 1))
            radius = max_dist + radius_extra
            return f"{cog[0]:.3f} {cog[1]:.3f} {cog[2]:.3f} {radius:.3f}"

        u = vg.Utils.create_mda_universe_quiet(path_pdb, path_traj)
        atoms = u.select_atoms(query)
        if len(atoms) == 0:
            raise ValueError(f"No atoms found for query: {query}")

        info = ' '.join(get_sphere_info() for _ in u.trajectory)
        vg.Utils.delete_traj_locks(path_traj)
        return info


# //////////////////////////////////////////////////////////////////////////////
