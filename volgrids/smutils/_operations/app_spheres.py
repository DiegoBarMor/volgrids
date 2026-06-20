import numpy as np
from pathlib import Path

import volgrids as vg
import volgrids.smiffer as sm
import volgrids.smutils as su
from volgrids._vendors import freyacli as fy

# //////////////////////////////////////////////////////////////////////////////
class AppSpheres(vg.AppSubcommand):
    # -------------------------------------------------------------------------- UI SECTION
    def run(self):
        operation = self.main.subcommands.pop(0)
        if operation == "find": return self.app_run_find()
        if operation == "grid": return self.app_run_grid()
        raise ValueError(f"Unknown operation: {operation}")


    # --------------------------------------------------------------------------
    def app_run_find(self):
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_traj = self.main.get_arg_path("path_traj",
            assertion = fy.PathAssertion.FILE_IN, allow_none = True
        )
        query = self.main.get_arg_str("query")
        radius_extra = self.main.get_arg_float("radius_extra")

        print(su.AppSpheres.find(path_in, path_traj, query, radius_extra))


    # --------------------------------------------------------------------------
    def app_run_grid(self):
        path_in  = self.main.get_arg_path("path_in", assertion = fy.PathAssertion.FILE_IN)
        path_traj = self.main.get_arg_path("path_traj",
            assertion = fy.PathAssertion.FILE_IN, allow_none = True
        )
        folder_out = self.main.get_arg_path("folder_out",
            assertion = fy.PathAssertion.DIR_OUT, default = path_in.parent
        )
        spheres_flat = self.main.get_arg_float("sphere", is_list = True)

        try: spheres = sm.SphereInfo.sphere_list(spheres_flat)
        except ValueError as e: self.main.help_and_exit(1, f"{e}")

        self.main.load_configs(vg, sm, su) # needed for loading sm.GRID_FORMAT_OUTPUT (used by Smif.save_data)

        su.AppSpheres.grid(path_in, path_traj, folder_out, spheres)
        vg.Utils.delete_traj_locks(path_traj)


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


    # --------------------------------------------------------------------------
    @staticmethod
    def grid(path_pdb: Path, path_traj: Path|None, folder_out: Path, spheres: list[sm.SphereInfo]) -> None:
        ms = sm.MolSystem(path_pdb, path_traj)
        nframes = ms.system.trajectory.n_frames

        sm.SphereInfo.assert_sphere_list(spheres, nframes)

        grid = vg.Grid(ms.box, dtype = bool)
        for i,sphere in enumerate(spheres):
            ms.switch_frame(i)
            timer = vg.Timer(f"...>>> Frame {ms.frame}/{nframes}")
            timer.start()

            grid.reset()
            kernel = vg.KernelSphere(sphere.radius, ms.get_deltas(), bool)
            kernel.stamp(grid, sphere.get_pos())
            path_out = folder_out / f"{ms.molname}.sphere.cmap"
            key_out = f"sphere.{i:04}"
            sm.Smif.save_data(grid, ms, path_out, key_out) # [WIP] check this

            timer.end()


# //////////////////////////////////////////////////////////////////////////////
