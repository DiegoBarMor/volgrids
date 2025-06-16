import numpy as np
import volgrids.vgrids as vg
import volgrids.smiffer as sm

# //////////////////////////////////////////////////////////////////////////////
class GridTrimmer(vg.Grid):
    def __init__(self, ms: "sm.SmifferMolecularSystem", trimming_dist):
        super().__init__(ms, dtype = bool)

        if sm.DO_TRIMMING_OCCUPANCY:
            self._trim_occupancy(trimming_dist)

        if sm.DO_TRIMMING_SPHERE and ms.meta.do_ps:
            self._trim_sphere()

        if sm.DO_TRIMMING_RNDS and ms.meta.do_ps:
            self._trim_rnds()

        if sm.SAVE_CACHED_MASK:
            self.save_data()


    # --------------------------------------------------------------------------
    def get_type(self):
        return "mask"


    # --------------------------------------------------------------------------
    def apply_trimming(self, vgrid: "vg.Grid"):
        vgrid.grid[self.grid] = 0


    # --------------------------------------------------------------------------
    def _trim_occupancy(self, radius):
        sk = vg.KernelSphere(radius, self.ms.deltas, bool)
        sk.link_to_grid(self.grid, self.ms.minCoords)
        for a in self.ms.relevant_atoms:
            sk.stamp(a.position)


    # --------------------------------------------------------------------------
    def _trim_sphere(self):
        coords = vg.get_coords_array(self.ms.resolution, self.ms.deltas, self.ms.minCoords)
        shifted_coords = coords - self.ms.cog
        dist_from_cog = vg.get_norm(shifted_coords)
        self.grid[dist_from_cog > self.ms.radius] = True


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
                if self.grid[neigh]: continue

                queue.add(neigh)

        self.grid[np.logical_not(visited)] = True


# //////////////////////////////////////////////////////////////////////////////
