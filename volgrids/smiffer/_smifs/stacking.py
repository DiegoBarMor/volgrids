import numpy as np
from collections import defaultdict

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifStacking(sm.Smif):
    def populate_grid(self, grid: vg.Grid) -> None:
        grid.reset()
        kernel = vg.KernelGaussianBivariateAngleDist(
            radius = sm.MU_DIST_STACKING + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_STACKING,
            deltas = self.ms.get_deltas(), dtype = vg.FLOAT_DTYPE, params = sm.PARAMS_STACK
        )
        for atoms_plane in self.iter_particles():
            cog, normal = self.get_cog_normal(atoms_plane)
            kernel.recalculate_kernel(normal, is_stacking = True)
            kernel.stamp(grid, cog, multiply_by = sm.ENERGY_SCALE)


    # --------------------------------------------------------------------------
    def iter_particles(self):
        import MDAnalysis as mda

        resname_to_ids = defaultdict(set)
        atoms = self.ms.get_relevant_atoms()

        if not len(atoms): return
        try:
            _ = atoms[0].chainID
            skip_chainID = False
        except mda.exceptions.NoDataError:
            skip_chainID = True


        for a in atoms:
            chain = None if skip_chainID else a.chainID
            resname_to_ids[a.resname.upper()].add((a.resid, chain))

        for resname,res_infos in resname_to_ids.items():
            lst_planes_atoms = self.ms.chemtable.get_names_stacking(resname)
            if lst_planes_atoms is None: continue

            for resid,chain in res_infos:
                for plane_atoms in lst_planes_atoms:
                    sel = f"resid {resid} and name {plane_atoms}"
                    if chain: sel += f" and chainID {chain}"

                    atoms_plane = atoms.select_atoms(sel)
                    if len(atoms_plane) < 3: continue # include rings even if they're not completely inside the grid's boundaries

                    yield atoms_plane


    # --------------------------------------------------------------------------
    @staticmethod
    def get_cog_normal(atoms_plane) -> tuple[np.ndarray, np.ndarray]:
        cog = atoms_plane.center_of_geometry()
        a,b,c = atoms_plane.positions[:3]
        u = vg.Math.normalize(b - a)
        v = vg.Math.normalize(c - a)
        normal = vg.Math.normalize(np.cross(u, v))
        return cog, normal


# //////////////////////////////////////////////////////////////////////////////
