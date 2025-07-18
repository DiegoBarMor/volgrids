import numpy as np
from collections import defaultdict

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class SmifStacking(sm.Smif):
    def populate_grid(self):
        radius = sm.MU_STACKING[1] + sm.GAUSSIAN_KERNEL_SIGMAS * sm.SIGMA_DIST_STACKING
        kernel = vg.KernelGaussianMultivariate(radius, self.ms.deltas, vg.FLOAT_DTYPE)

        kernel.link_to_grid(self.grid, self.ms.minCoords)
        for res_atoms in self.iter_particles():
            cog = res_atoms.center_of_geometry()
            a,b,c = res_atoms.positions[:3]
            u = vg.Math.normalize(b - a)
            v = vg.Math.normalize(c - a)
            normal = vg.Math.normalize(np.cross(u, v))

            kernel.recalculate_kernel(normal, sm.MU_STACKING, sm.COV_INV_STACKING, isStacking = True)
            kernel.stamp(cog, multiplication_factor = sm.ENERGY_SCALE)


    # --------------------------------------------------------------------------
    def iter_particles(self):
        resname_to_ids = defaultdict(set)
        atoms = self.ms.get_relevant_atoms()

        for a in atoms:
            resname_to_ids[a.resname.upper()].add((a.resid, a.chainID))

        for resname,res_infos in resname_to_ids.items():
            aromatic_atoms = self.ms.chemtable.get_names_stacking(resname)
            if aromatic_atoms is None: continue

            for resid,chain in res_infos:
                sel = f"resid {resid} and name {aromatic_atoms}"
                if chain: sel += f" and chainID {chain}"
                res_atoms = atoms.select_atoms(sel)
                if len(res_atoms) >= 3: # include rings even if they're not completely inside the PS
                    yield res_atoms


# //////////////////////////////////////////////////////////////////////////////
