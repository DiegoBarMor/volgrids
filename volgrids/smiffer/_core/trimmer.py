import numpy as np

import volgrids as vg
import volgrids.smiffer as smf

# //////////////////////////////////////////////////////////////////////////////
class Trimmer:
    def __init__(self, mm: "smf.MoleculeManager"):
        self.mm: "smf.MoleculeManager" = mm
        self.cavfinder: smf.CavityFinder = smf.CavityFinder()

        self._mask_common: vg.Grid = None
        self._mask_specific: vg.Grid = None
        self._current_key: str = ""
        self._distances = {
            "small": vg.CFG.trim_occ_dist_short,
            "mid"  : vg.CFG.trim_occ_dist_mid,
            "large": vg.CFG.trim_occ_dist_long,
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
        return vg.CFG.smif_hphil


    # --------------------------------------------------------------------------
    @classmethod
    def should_do_trim_mid(cls) -> bool:
        return any((
            vg.CFG.smif_stk, vg.CFG.smif_hba, vg.CFG.smif_hbd, vg.CFG.smif_hphob,
            vg.CFG.trim_save, cls.should_do_cavities()
        ))


    # --------------------------------------------------------------------------
    @classmethod
    def should_do_trim_large(cls) -> bool:
        return vg.CFG.smif_apbs


    # --------------------------------------------------------------------------
    @staticmethod
    def should_do_cavities() -> bool:
        return any((
            vg.CFG.trim_cavities, vg.CFG.cav_save,
            vg.CFG.cav_weight != 0.0
        )) and vg.CFG.trim_occupancy


    # --------------------------------------------------------------------------
    def _should_run_mask_common(self) -> bool:
        return self.mm.do_use_sphere and self._mask_common is None


    # --------------------------------------------------------------------------
    def _should_run_mask_specific(self, key: str) -> bool:
        if key == self._current_key: return False
        self._current_key = key
        return vg.CFG.trim_occupancy


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
        self._mask_common = vg.Grid(self.mm.box, dtype = bool)
        if vg.CFG.trim_faraway: self._trim_faraway()
        if vg.CFG.trim_sphere: self._trim_sphere()
        if vg.CFG.trim_rnds: self._trim_rnds()


    # --------------------------------------------------------------------------
    def _run_specific(self, dist: float):
        if self._mask_specific is None:
            self._mask_specific = vg.Grid(self.mm.box, dtype = bool)
        else:
            self._mask_specific.reset()

        self._trim_occupancies(dist)


    # --------------------------------------------------------------------------
    def _run_cavities(self):
        """must be called immediately after `_trim_occupancies`, before any other trimming operations"""

        self.cavfinder.populate_cavities_grid(self._mask_specific)

        if not vg.CFG.trim_cavities: return

        if self._mask_common is None:
            self._mask_common = vg.Grid(self.mm.box, dtype = bool)

        self._mask_common.arr |= (self.cavfinder.grid.arr < vg.CFG.cav_threshold)


    # --------------------------------------------------------------------------
    def _trim_occupancies(self, radius: float):
        kernel = vg.KernelSphere(radius, self.mm.get_deltas(), bool)
        for a in self.mm.get_all_queried_atoms(use_custom = False):
            kernel.stamp(self._mask_specific, a.get_position_numpy())


    # --------------------------------------------------------------------------
    def _trim_sphere(self):
        coords = vg.Math.get_coords_array(self.mm.get_resolution(), self.mm.get_deltas(), self.mm.get_min_coords())
        shifted_coords = coords - self.mm.get_cog()
        dist_from_cog = vg.Math.get_norm(shifted_coords)
        self._mask_common.arr[dist_from_cog > self.mm.get_radius()] = True


    # --------------------------------------------------------------------------
    def _trim_rnds(self):
        """
        Perform a random search to remove isolated regions.
        Can be problematic (e.g. slow, aggressive trimming); use with caution.
        """
        visited = np.zeros(self.mm.get_resolution(), dtype = bool)

        directions = np.array([[i,j,k] for i in range(-1,2) for j in range(-1,2) for k in range(-1,2) if i&j&k])

        xres, yres, zres = self.mm.get_resolution()
        xcog, ycog, zcog = np.floor(self.mm.get_resolution() / 2).astype(int)
        cog_cube = set((x,y,z)
            for x in range(xcog - vg.CFG.trim_rnds_cube_radius, xcog + vg.CFG.trim_rnds_cube_radius + 1)
            for y in range(ycog - vg.CFG.trim_rnds_cube_radius, ycog + vg.CFG.trim_rnds_cube_radius + 1)
            for z in range(zcog - vg.CFG.trim_rnds_cube_radius, zcog + vg.CFG.trim_rnds_cube_radius + 1)
        )
        queue = cog_cube.copy()

        search_dist = np.full(self.mm.get_resolution(), np.inf)
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
                if search_dist[neigh] > vg.CFG.trim_rnds_max_dist: continue
                if visited[neigh]: continue
                if self._mask_common.arr[neigh]: continue

                queue.add(neigh)

        self._mask_common.arr[np.logical_not(visited)] = True


    # --------------------------------------------------------------------------
    def _trim_faraway(self):
        grid = self._mask_common.copy()
        kernel = vg.KernelSphere(vg.CFG.trim_faraway_dist, self.mm.get_deltas(), bool)
        for a in self.mm.get_all_queried_atoms(use_custom = False):
            kernel.stamp(grid, a.get_position_numpy())
        self._mask_common.arr[~grid.arr] = True


# //////////////////////////////////////////////////////////////////////////////
