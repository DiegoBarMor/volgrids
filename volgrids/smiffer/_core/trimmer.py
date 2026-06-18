import numpy as np

import volgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class Trimmer:
    def __init__(self, ms: "sm.MolSystem"):
        self.ms: "sm.MolSystem" = ms
        self.cavfinder: sm.CavityFinder = sm.CavityFinder()

        self._mask_common: vg.Grid = None
        self._mask_specific: vg.Grid = None
        self._current_key: str = ""
        self._distances = {
            "small": sm.TRIMMING_DIST_SMALL,
            "mid"  : sm.TRIMMING_DIST_MID,
            "large": sm.TRIMMING_DIST_LARGE,
        }


    # --------------------------------------------------------------------------
    def trim(self, grid: "vg.Grid", key: str):
        """Removes grid points wherever the mask (for the given `key`) is `True`."""

        if self._should_run_mask_common(): # will only run in the first call to `trim`
            self._run_common()

        if self._should_run_mask_specific(key): # can run in multiple calls to `trim`
            self._run_specific(dist = self._distances[key])

        if self._should_run_cavities(): # will only run in the first call to `trim`
            self._run_cavities()

        mask = self._get_mask_merged()
        if mask is None: return

        grid.arr[mask.arr] = 0
        self.cavfinder.apply_cavities_weighting(grid)


    # --------------------------------------------------------------------------
    def get_mask(self) -> vg.Grid|None:
        """Returns the merged mask (if any) from the trimming operations."""
        return self._get_mask_merged()


    # --------------------------------------------------------------------------
    def run_for_saving(self, key: str) -> None:
        """
        Trimming operations could be skipped if the only intention is to save the trimming/cavities grid.
        These method ensures that the necessary trimming operations are performed in this case.
        """
        if self._should_run_mask_specific(key): self._run_specific(dist = self._distances[key])
        if self._should_run_cavities(): self._run_cavities()


    # --------------------------------------------------------------------------
    @classmethod
    def should_do_trim_small(cls) -> bool:
        return sm.DO_SMIF_HYDROPHILIC


    # --------------------------------------------------------------------------
    @classmethod
    def should_do_trim_mid(cls) -> bool:
        return any((
            sm.DO_SMIF_STACKING, sm.DO_SMIF_HBA, sm.DO_SMIF_HBD, sm.DO_SMIF_HYDROPHOBIC,
            sm.SAVE_TRIMMING_MASK, cls.should_do_cavities()
        ))


    # --------------------------------------------------------------------------
    @classmethod
    def should_do_trim_large(cls) -> bool:
        return sm.DO_SMIF_APBS


    # --------------------------------------------------------------------------
    @staticmethod
    def should_do_cavities() -> bool:
        return any((
            sm.DO_TRIMMING_CAVITIES, sm.SAVE_CAVITIES,
            sm.CAVITIES_WEIGHT != 0.0
        )) and sm.DO_TRIMMING_OCCUPANCY


    # --------------------------------------------------------------------------
    def _should_run_mask_common(self) -> bool:
        return self.ms.do_ps and self._mask_common is None


    # --------------------------------------------------------------------------
    def _should_run_mask_specific(self, key: str) -> bool:
        if key == self._current_key: return False
        self._current_key = key
        return sm.DO_TRIMMING_OCCUPANCY


    # --------------------------------------------------------------------------
    def _should_run_cavities(self) -> bool:
        return self.should_do_cavities() and not self.cavfinder.has_data()


    # --------------------------------------------------------------------------
    def _get_mask_merged(self) -> vg.Grid|None:
        if (self._mask_common is None) and (self._mask_specific is None): return
        if self._mask_specific is None: return self._mask_common
        if self._mask_common is None: return self._mask_specific

        self._mask_specific.arr |= self._mask_common.arr
        return self._mask_specific


    # --------------------------------------------------------------------------
    def _run_common(self):
        self._mask_common = vg.Grid(self.ms.box, dtype = bool)
        if sm.DO_TRIMMING_FARAWAY: self._trim_faraway()
        if sm.DO_TRIMMING_SPHERE: self._trim_sphere()
        if sm.DO_TRIMMING_RNDS: self._trim_rnds()


    # --------------------------------------------------------------------------
    def _run_specific(self, dist: float):
        if self._mask_specific is None:
            self._mask_specific = vg.Grid(self.ms.box, dtype = bool)
        else:
            self._mask_specific.reset()

        self._trim_occupancies(dist)


    # --------------------------------------------------------------------------
    def _run_cavities(self):
        """must be called immediately after `_trim_occupancies`, before any other trimming operations"""

        self.cavfinder.populate_cavities_grid(self._mask_specific)

        if not sm.DO_TRIMMING_CAVITIES: return

        if self._mask_common is None:
            self._mask_common = vg.Grid(self.ms.box, dtype = bool)

        self._mask_common.arr |= (self.cavfinder.grid.arr < sm.TRIMMING_CAVITIES_THRESHOLD)


    # --------------------------------------------------------------------------
    def _trim_occupancies(self, radius: float):
        kernel = vg.KernelSphere(radius, self.ms.get_deltas(), bool)
        for a in self.ms.get_all_atoms(use_custom = False):
            kernel.stamp(self._mask_specific, a.position)


    # --------------------------------------------------------------------------
    def _trim_sphere(self):
        coords = vg.Math.get_coords_array(self.ms.get_resolution(), self.ms.get_deltas(), self.ms.get_min_coords())
        shifted_coords = coords - self.ms.get_cog()
        dist_from_cog = vg.Math.get_norm(shifted_coords)
        self._mask_common.arr[dist_from_cog > self.ms.get_radius()] = True


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
                if self._mask_common.arr[neigh]: continue

                queue.add(neigh)

        self._mask_common.arr[np.logical_not(visited)] = True


    # --------------------------------------------------------------------------
    def _trim_faraway(self):
        grid = self._mask_common.copy()
        kernel = vg.KernelSphere(sm.TRIM_FARAWAY_DIST, self.ms.get_deltas(), bool)
        for a in self.ms.get_all_atoms(use_custom = False):
            kernel.stamp(grid, a.position)
        self._mask_common.arr[~grid.arr] = True


# //////////////////////////////////////////////////////////////////////////////
