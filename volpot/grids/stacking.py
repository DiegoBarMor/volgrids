import volpot
import numpy as np
from collections import defaultdict

from settings import MU_STACKING, GAUSSIAN_KERNEL_SIGMAS, SIGMA_DIST_STACKING, COV_INV_STACKING


# //////////////////////////////////////////////////////////////////////////////
class SPG_Stacking(volpot.StatisticalPotentialGrid):
    POTENTIAL_TYPE = "stacking"

    def populate_grid(self):
        radius = MU_STACKING[1] + GAUSSIAN_KERNEL_SIGMAS * SIGMA_DIST_STACKING
        gk = volpot.MultiGaussianKernel(radius, self.ms.deltas, np.float32)

        gk.link_to_grid(self.grid, self.ms.minCoords)
        for res_atoms in self.iter_particles():
            cog = res_atoms.center_of_geometry()
            a,b,c = res_atoms.positions[:3]
            u = volpot.normalize(b - a)
            v = volpot.normalize(c - a)
            normal = volpot.normalize(np.cross(u, v))

            gk.recalculate_kernel(normal, MU_STACKING, COV_INV_STACKING, isStacking = True)
            gk.stamp(cog)


    def iter_particles(self):
        resname_to_ids = defaultdict(set)
        aromatic_dict = volpot.aromatic_bases if self.ms.isNucleic else volpot.aromatic_aminos

        for a in self.ms.relevant_atoms:
            resname_to_ids[a.resname].add((a.resid, a.chainID))

        for resname,res_infos in resname_to_ids.items():
            aromatic_ring = aromatic_dict.get(resname)
            if aromatic_ring is None: continue
            aromatic_atoms, ring_size = aromatic_ring

            for resid,chain in res_infos:
                sel = f"resid {resid} and name {aromatic_atoms}"
                if chain: sel += f" and chainID {chain}"
                res_atoms = self.ms.relevant_atoms.select_atoms(sel)
                if len(res_atoms) >= 3: # include rings even if they're not completely inside the PS
                    yield res_atoms


# //////////////////////////////////////////////////////////////////////////////
