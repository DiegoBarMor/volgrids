import numpy as np

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Trimmer:
    KEY_INIT_COMMON_MASK = "mid" # the common mask is initialized by copying this specific mask

    def __init__(self, ms: "sm.MolSystemSmiffer", **distances):
        self.ms: "sm.MolSystemSmiffer" = ms

        self.distances: dict[str, float] = distances
        self.common_mask: vg.Grid = None
        self.specific_masks: dict[str, vg.Grid] = None


    # --------------------------------------------------------------------------
    @classmethod
    def init_infer_dists(cls, ms: "sm.MolSystemSmiffer") -> "sm.Trimmer":
        trimming_dists = {}
        if any((
            sm.DO_SMIF_PP, sm.DO_SMIF_HBA_PP, sm.DO_SMIF_HBD_PP,
            sm.DO_SMIF_STACKING_PP, sm.DO_SMIF_HYDRO_PP
        )):
            trimming_dists["tiny"] = sm.TRIMMING_DIST_TINY

        if sm.DO_SMIF_HYDROPHILIC:
            trimming_dists["small"] = sm.TRIMMING_DIST_SMALL

        if any((
            sm.DO_SMIF_STACKING, sm.DO_SMIF_HBA, sm.DO_SMIF_HBD, sm.DO_SMIF_HYDROPHOBIC,
            sm.SAVE_TRIMMING_MASK, cls._should_do_cavities()
        )):
            trimming_dists["mid"] = sm.TRIMMING_DIST_MID

        if sm.DO_SMIF_APBS:
            trimming_dists["large"] = sm.TRIMMING_DIST_LARGE

        return cls(ms, **trimming_dists)


    # --------------------------------------------------------------------------
    def reset_masks(self):
        self.common_mask: vg.Grid = None
        self.specific_masks = {k : vg.Grid(self.ms.box, dtype = bool) for k in self.distances.keys()}


    # --------------------------------------------------------------------------
    def trim(self, cavfinder: "sm.CavityFinder"):
        self.reset_masks()

        if sm.DO_TRIMMING_OCCUPANCY:
            self._trim_occupancies()
            self._trim_cavities(cavfinder)

        if self._should_use_common_mask():
            self._init_common_mask()
            self._run_common_mask_operations()
            self._apply_common_mask_to_specific_masks()
            self._discard_common_mask()


    # --------------------------------------------------------------------------
    def mask_grid(self, smif: "sm.Smif", key: str):
        """Removes `smif` points wherever the mask (for the given `key`) is `True`."""
        smif.grid.arr[self.get_mask(key).arr] = 0


    # --------------------------------------------------------------------------
    def get_mask(self, key: str) -> vg.Grid:
        return self.specific_masks[key]


    # --------------------------------------------------------------------------
    def _get_smallest_mask(self) -> vg.Grid:
        key_smallest = min(self.specific_masks.keys(), key = lambda k: self.distances[k])
        return self.specific_masks[key_smallest]


    # --------------------------------------------------------------------------
    @staticmethod
    def _should_do_cavities() -> bool:
        return any((
            sm.DO_TRIMMING_CAVITIES, sm.SAVE_CAVITIES,
            sm.CAVITIES_WEIGHT != 0.0
        ))

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
            mask.arr |= self.common_mask.arr


    # --------------------------------------------------------------------------
    def _discard_common_mask(self):
        del self.common_mask
        self.common_mask = None


    # --------------------------------------------------------------------------
    def _trim_occupancies(self):
        for k,radius in self.distances.items():
            mask = self.specific_masks[k]
            kernel = vg.KernelSphere(radius, self.ms.get_deltas(), bool)
            for a in self.ms.get_relevant_atoms(use_custom = False, extra_dist = radius):
                kernel.stamp(mask, a.position)


    # --------------------------------------------------------------------------
    def _trim_cavities(self, cavfinder: "sm.CavityFinder"):
        """must be called immediately after `_trim_occupancies`, before any other trimming operations"""

        if not self._should_do_cavities(): return

        cavfinder.populate_cavities_grid(self._get_smallest_mask())

        if sm.DO_TRIMMING_CAVITIES:
            for mask in self.specific_masks.values():
                mask.arr |= (cavfinder.grid.arr < sm.TRIMMING_CAVITIES_THRESHOLD)


    # --------------------------------------------------------------------------
    def _trim_sphere(self):
        coords = vg.Math.get_coords_array(self.ms.get_resolution(), self.ms.get_deltas(), self.ms.get_min_coords())
        shifted_coords = coords - self.ms.get_cog()
        dist_from_cog = vg.Math.get_norm(shifted_coords)
        self.common_mask.arr[dist_from_cog > self.ms.get_radius()] = True


    # --------------------------------------------------------------------------
    def _trim_rnds(self):
        """
        Perform a random search to remove isolated regions.
        Can be problematic (e.g. slow, aggressive trimming); use with caution.
        """
        visited = np.zeros(self.ms.get_resolution(), dtype = bool)

        directions = np.array([[i,j,k] for i in range(-1,2) for j in range(-1,2) for k in range(-1,2) if i&j&k])

        xres, yres, zres = self.ms.get_resolution()
        xcog, ycog, zcog = np.floor(self.ms.get_resolution() / 2).astype(int)
        cog_cube = set((x,y,z)
            for x in range(xcog - sm.COG_CUBE_RADIUS, xcog + sm.COG_CUBE_RADIUS + 1)
            for y in range(ycog - sm.COG_CUBE_RADIUS, ycog + sm.COG_CUBE_RADIUS + 1)
            for z in range(zcog - sm.COG_CUBE_RADIUS, zcog + sm.COG_CUBE_RADIUS + 1)
        )
        queue = cog_cube.copy()

        search_dist = np.full(self.ms.get_resolution(), np.inf)
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
                if self.common_mask.arr[neigh]: continue

                queue.add(neigh)

        self.common_mask.arr[np.logical_not(visited)] = True


    # --------------------------------------------------------------------------
    def _trim_faraway(self):
        grid = self.common_mask.copy()
        kernel = vg.KernelSphere(sm.TRIM_FARAWAY_DIST, self.ms.get_deltas(), bool)
        for a in self.ms.get_relevant_atoms(use_custom = False, extra_dist = sm.TRIM_FARAWAY_DIST):
            kernel.stamp(grid, a.position)
        self.common_mask.arr[~grid.arr] = True


# //////////////////////////////////////////////////////////////////////////////
