import numpy as np

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Trimmer:
    KEY_INIT_COMMON_MASK = "mid" # the common mask is initialized by copying this specific mask

    def __init__(self, ms: "sm.MolSystemSmiffer", **distances):
        self.ms: "sm.MolSystemSmiffer" = ms

        self.distances = distances
        self.common_mask: vg.Grid = None
        self.specific_masks = {k : vg.Grid(ms, dtype = bool) for k in distances.keys()}


    # --------------------------------------------------------------------------
    @classmethod
    def init_infer_dists(cls, ms: "sm.MolSystemSmiffer") -> "sm.Trimmer":
        trimming_dists = {}
        if sm.DO_SMIF_HYDROPHILIC:
            trimming_dists["small"] = sm.TRIMMING_DIST_SMALL

        if (
            sm.DO_SMIF_STACKING or
            sm.DO_SMIF_HBA or sm.DO_SMIF_HBD or
            sm.DO_SMIF_HYDROPHOBIC or sm.SAVE_TRIMMING_MASK
        ):
            trimming_dists["mid"] = sm.TRIMMING_DIST_MID

        if sm.DO_SMIF_APBS:
            trimming_dists["large"] = sm.TRIMMING_DIST_LARGE

        return cls(ms, **trimming_dists)


    # --------------------------------------------------------------------------
    def trim(self):
        if sm.DO_TRIMMING_OCCUPANCY:
            self._trim_occupancies()

        if self._should_use_common_mask():
            self._init_common_mask()
            self._run_common_mask_operations()
            self._apply_common_mask_to_specific_masks()
            self._discard_common_mask()


    # --------------------------------------------------------------------------
    def mask_grid(self, smif: "vg.Grid", key: str):
        smif.grid[self.get_mask(key).grid] = 0


    # --------------------------------------------------------------------------
    def get_mask(self, key: str) -> vg.Grid:
        return self.specific_masks[key]


    # --------------------------------------------------------------------------
    def _should_use_common_mask(self) -> bool:
        return self.ms.do_ps


    # --------------------------------------------------------------------------
    def _init_common_mask(self):
        self.common_mask = self.specific_masks[self.KEY_INIT_COMMON_MASK].copy()


    # --------------------------------------------------------------------------
    def _run_common_mask_operations(self):
        if sm.DO_TRIMMING_FARAWAY:
            self._trim_faraway()

        if sm.DO_TRIMMING_SPHERE:
            self._trim_sphere()

        if sm.DO_TRIMMING_RNDS:
            self._trim_rnds()


    # --------------------------------------------------------------------------
    def _apply_common_mask_to_specific_masks(self):
        for mask in self.specific_masks.values():
            mask.grid |= self.common_mask.grid


    # --------------------------------------------------------------------------
    def _discard_common_mask(self):
        del self.common_mask
        self.common_mask = None


    # --------------------------------------------------------------------------
    def _trim_occupancies(self):
        for k,radius in self.distances.items():
            mask = self.specific_masks[k]
            kernel = vg.KernelSphere(radius, self.ms.deltas, bool)
            kernel.link_to_grid(mask.grid, self.ms.minCoords)

            for a in self.ms.get_relevant_atoms_broad(radius):
                kernel.stamp(a.position)


    # --------------------------------------------------------------------------
    def _trim_sphere(self):
        coords = vg.Math.get_coords_array(self.ms.resolution, self.ms.deltas, self.ms.minCoords)
        shifted_coords = coords - self.ms.cog
        dist_from_cog = vg.Math.get_norm(shifted_coords)
        self.common_mask.grid[dist_from_cog > self.ms.radius] = True


    # --------------------------------------------------------------------------
    def _trim_rnds(self):
        """perform a random search to remove isolated regions"""
        visited = np.zeros(self.ms.resolution, dtype = bool)

        directions = np.array([[i,j,k] for i in range(-1,2) for j in range(-1,2) for k in range(-1,2) if i&j&k])

        xres, yres, zres = self.ms.resolution
        xcog, ycog, zcog = np.floor(self.ms.resolution / 2).astype(int)
        cog_cube = set((x,y,z)
            for x in range(xcog - sm.COG_CUBE_RADIUS, xcog + sm.COG_CUBE_RADIUS + 1)
            for y in range(ycog - sm.COG_CUBE_RADIUS, ycog + sm.COG_CUBE_RADIUS + 1)
            for z in range(zcog - sm.COG_CUBE_RADIUS, zcog + sm.COG_CUBE_RADIUS + 1)
        )
        queue = cog_cube.copy()

        search_dist = np.full(self.ms.resolution, np.inf)
        for point in cog_cube: search_dist[point] = 0

        while queue:
            ### "random search" because popping from a set can be unpredictable
            i,j,k = node = queue.pop()
            visited[node] = True

            for dx,dy,dz in directions:
                ni = i + dx
                if not (0 <= ni < xres): continue

                nj = j + dy
                if not (0 <= nj < yres): continue

                nk = k + dz
                if not (0 <= nk < zres): continue

                neigh = ni,nj,nk
                search_dist[neigh] = min(search_dist[node] + 1, search_dist[neigh])
                if search_dist[neigh] > sm.MAX_RNDS_DIST: continue
                if visited[neigh]: continue
                if self.common_mask.grid[neigh]: continue

                queue.add(neigh)

        self.common_mask.grid[np.logical_not(visited)] = True


    # --------------------------------------------------------------------------
    def _trim_faraway(self):
        arr = np.zeros_like(self.common_mask.grid, dtype = bool)
        kernel = vg.KernelSphere(sm.TRIM_FARAWAY_DIST, self.ms.deltas, bool)
        kernel.link_to_grid(arr, self.ms.minCoords)
        for a in self.ms.get_relevant_atoms_broad(sm.TRIM_FARAWAY_DIST):
            kernel.stamp(a.position)

        self.common_mask.grid[~arr] = True


# //////////////////////////////////////////////////////////////////////////////
