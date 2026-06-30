import numpy as np

import volgrids as vg
import volgrids.smiffer as smf
from volgrids._vendors import molsimple as ms

# //////////////////////////////////////////////////////////////////////////////
class SmifStacking(smf.Smif):
    def populate_grid(self, grid: vg.Grid) -> None:
        grid.reset()
        kernel = vg.KernelGaussianBivariateAngleDist(
            radius = vg.CFG.param_stk_dist_mu + vg.CFG.misc_kernel_gaussian_sigmas * smf.SIGMA_DIST_STACKING,
            deltas = self.mm.get_deltas(), dtype = vg.FLOAT_DTYPE, params = smf.PARAMS_STACK
        )
        for atoms_plane in self.iter_particles():
            cog, normal = self.get_cog_normal(atoms_plane)
            kernel.recalculate_kernel(normal, is_stacking = True)
            kernel.stamp(grid, cog, multiply_by = vg.CFG.param_stk_scale)


    # --------------------------------------------------------------------------
    def iter_particles(self):
        atoms = self.mm.get_relevant_queried_atoms()
        if not len(atoms): return

        residues = atoms.split_residues()

        for residue in residues:
            resname = residue[0].resname.upper()
            lst_names_planes = self.mm.chemtable.get_names_stacking(resname)

            if lst_names_planes is None: continue
            for names_plane in lst_names_planes:
                atoms_plane = residue.select_name(*names_plane.split())
                if len(atoms_plane) < 3: continue # include rings even if they're not completely inside the grid's boundaries

                yield atoms_plane


    # --------------------------------------------------------------------------
    @staticmethod
    def get_cog_normal(atoms_plane: ms.ParticleGroup) -> tuple[np.ndarray, np.ndarray]:
        coords = atoms_plane.get_positions_numpy()
        cog = np.mean(coords, axis = 0)
        a,b,c = coords[:3]
        u = vg.Math.normalize(b - a)
        v = vg.Math.normalize(c - a)
        normal = vg.Math.normalize(np.cross(u, v))
        return cog, normal


# //////////////////////////////////////////////////////////////////////////////
