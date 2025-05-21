import numpy as np
import volpot as vp
from abc import ABC, abstractmethod

# //////////////////////////////////////////////////////////////////////////////
class GridTrimmer(ABC, vp.Grid):
    def __init__(self, ms, trimming_dist):
        super().__init__(ms, dtype = bool)

        if vp.DO_TRIMMING_SPHERE:    self._trim_sphere()
        if vp.DO_TRIMMING_OCCUPANCY: self._trim_occupancy(trimming_dist)
        if vp.DO_TRIMMING_RNDS:      self._trim_rnds()

        if vp.SAVE_CACHED_MASK:
            self.pack_data()
            self.save_data()

    def get_type(self):
        return "mask"

    def apply_trimming(self, potential_grid):
        potential_grid.grid[self.grid] = 0


    def _trim_occupancy(self, radius):
        sk = vp.KernelSphere(radius, self.ms.deltas, bool)
        sk.link_to_grid(self.grid, self.ms.minCoords)
        for a in self.ms.get_relevant_atoms_broad(radius):
            sk.stamp(a.position)

    @abstractmethod
    def _trim_sphere(self):
        return

    @abstractmethod
    def _trim_rnds(self):
        return


# //////////////////////////////////////////////////////////////////////////////
class TrimmerPocketSphere(GridTrimmer):
    def _trim_sphere(self):
        coords = vp.get_coords_array(self.ms.resolution, self.ms.deltas, self.ms.minCoords)
        shifted_coords = coords - self.ms.cog
        dist_from_cog = vp.get_norm(shifted_coords)
        self.grid[dist_from_cog > self.ms.radius] = True


    def _trim_rnds(self):
        """perform a random search to remove isolated regions"""
        visited = np.zeros(self.ms.resolution, dtype = bool)

        directions = np.array([[i,j,k] for i in range(-1,2) for j in range(-1,2) for k in range(-1,2) if i&j&k])

        xres, yres, zres = self.ms.resolution
        xcog, ycog, zcog = np.floor(self.ms.resolution / 2).astype(int)
        cog_cube = set((x,y,z)
            for x in range(xcog - vp.COG_CUBE_RADIUS, xcog + vp.COG_CUBE_RADIUS + 1)
            for y in range(ycog - vp.COG_CUBE_RADIUS, ycog + vp.COG_CUBE_RADIUS + 1)
            for z in range(zcog - vp.COG_CUBE_RADIUS, zcog + vp.COG_CUBE_RADIUS + 1)
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
                if search_dist[neigh] > vp.MAX_RNDS_DIST: continue
                if visited[neigh]: continue
                if self.grid[neigh]: continue

                queue.add(neigh)

        self.grid[np.logical_not(visited)] = True


# //////////////////////////////////////////////////////////////////////////////
class TrimmerWhole(GridTrimmer):
    def _trim_sphere(self):
        return

    def _trim_rnds(self):
        return


# //////////////////////////////////////////////////////////////////////////////
