import numpy as np
import volgrids as vg
from collections import defaultdict

# //////////////////////////////////////////////////////////////////////////////
class GridStacking(vg.GridSMIF):
    def get_type(self):
        return "stacking"

    def populate_grid(self):
        radius = vg.MU_STACKING[1] + vg.GAUSSIAN_KERNEL_SIGMAS * vg.SIGMA_DIST_STACKING
        gk = vg.KernelGaussianMultivariate(radius, self.ms.deltas, vg.FLOAT_DTYPE)

        gk.link_to_grid(self.grid, self.ms.minCoords)
        for res_atoms in self.iter_particles():
            cog = res_atoms.center_of_geometry()
            a,b,c = res_atoms.positions[:3]
            u = vg.normalize(b - a)
            v = vg.normalize(c - a)
            normal = vg.normalize(np.cross(u, v))

            gk.recalculate_kernel(normal, vg.MU_STACKING, vg.COV_INV_STACKING, isStacking = True)
            gk.stamp(cog)


    def iter_particles(self):
        resname_to_ids = defaultdict(set)
        aromatic_dict = vg.planar_rna if self.ms.isNucleic else vg.planar_prot
        atoms = self.ms.get_relevant_atoms()

        for a in atoms:
            resname_to_ids[a.resname.upper()].add((a.resid, a.chainID))

        for resname,res_infos in resname_to_ids.items():
            aromatic_atoms = aromatic_dict.get(resname)
            if aromatic_atoms is None: continue

            for resid,chain in res_infos:
                sel = f"resid {resid} and name {aromatic_atoms}"
                if chain: sel += f" and chainID {chain}"
                res_atoms = atoms.select_atoms(sel)
                if len(res_atoms) >= 3: # include rings even if they're not completely inside the PS
                    yield res_atoms


# //////////////////////////////////////////////////////////////////////////////
