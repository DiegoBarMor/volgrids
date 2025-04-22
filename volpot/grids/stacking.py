import numpy as np
import volpot as vp
from collections import defaultdict

# //////////////////////////////////////////////////////////////////////////////
class GridStacking(vp.GridSMIF):
    def get_type(self):
        return "stacking"

    def populate_grid(self):
        radius = vp.MU_STACKING[1] + vp.GAUSSIAN_KERNEL_SIGMAS * vp.SIGMA_DIST_STACKING
        gk = vp.KernelGaussianMultivariate(radius, self.ms.deltas, np.float32)

        gk.link_to_grid(self.grid, self.ms.minCoords)
        for res_atoms in self.iter_particles():
            cog = res_atoms.center_of_geometry()
            a,b,c = res_atoms.positions[:3]
            u = vp.normalize(b - a)
            v = vp.normalize(c - a)
            normal = vp.normalize(np.cross(u, v))

            gk.recalculate_kernel(normal, vp.MU_STACKING, vp.COV_INV_STACKING, isStacking = True)
            gk.stamp(cog)


    def iter_particles(self):
        resname_to_ids = defaultdict(set)
        aromatic_dict = vp.aromatic_bases if self.ms.isNucleic else vp.aromatic_aminos
        atoms = self.ms.get_relevant_atoms()

        for a in atoms:
            resname_to_ids[a.resname].add((a.resid, a.chainID))

        for resname,res_infos in resname_to_ids.items():
            aromatic_ring = aromatic_dict.get(resname)
            if aromatic_ring is None: continue
            aromatic_atoms, ring_size = aromatic_ring

            for resid,chain in res_infos:
                sel = f"resid {resid} and name {aromatic_atoms}"
                if chain: sel += f" and chainID {chain}"
                res_atoms = atoms.select_atoms(sel)
                if len(res_atoms) >= 3: # include rings even if they're not completely inside the PS
                    yield res_atoms


# //////////////////////////////////////////////////////////////////////////////
