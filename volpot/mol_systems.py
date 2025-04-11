import volpot
import numpy as np
from pathlib import Path
import MDAnalysis as mda
from gridData import Grid
from gridData import mrc

from settings import WARNING_GRID_SIZE, SAVE_CACHED_MASK, \
    DO_TRIMMING_SPHERE, DO_TRIMMING_OCCUPANCY, DO_TRIMMING_RNDS, \
    GRID_XRES_PS, GRID_YRES_PS, GRID_ZRES_PS, GRID_DX_PS, GRID_DY_PS, GRID_DZ_PS, \
    GRID_XRES_WM, GRID_YRES_WM, GRID_ZRES_WM, GRID_DX_WM, GRID_DY_WM, GRID_DZ_WM, \
    COG_CUBE_RADIUS, MAX_RNDS_DIST, EXTRA_BOX_SIZE, VERBOSE, DO_DX_OUTPUT, DO_MRC_OUTPUT


# //////////////////////////////////////////////////////////////////////////////
class MolecularSystem:
    def __init__(self, metadata, trimming_dist):
        self.resolution = np.zeros(3, dtype = int)
        self.minCoords = np.zeros(3)
        self.deltas = np.zeros(3)

        self.dx_comments = '\n'.join((
            "# OpenDX density file written by dx_io.DXIO.write_dx()",
            "# File format: http://opendx.sdsc.edu/docs/html/pages/usrgu068.htm#HDREDF",
            "# Data are embedded in the header and tied to the grid positions.",
            "# Data is written in C array order: In grid[x,y,z] the axis z is fastest",
            "# varying, then y, then finally x, i.e. z is the innermost loop.",
            '',
        ))
        self.dx_footer = '\n'.join((
            '',
            'attribute "dep" string "positions"',
            'object "density" class field',
            'component "positions" value 1',
            'component "connections" value 2',
            'component "data" value 3',
        ))

        self.path_pdb = Path(metadata["pdb"])
        self.path_apbs = Path(metadata["apbs"])
        self.path_metadata = Path(metadata["meta"])
        self.folder_potentials = Path(metadata["out"])

        self.name = self.path_pdb.stem
        self.isNucleic = bool(metadata["rna"])
        self.macro_query = f"{'nucleic' if self.isNucleic else 'protein'} and not (name H*)"

        self.system = mda.Universe(metadata["pdb"])
        self.metadata = metadata.copy()
        self.trimming_dist = trimming_dist

        ###############################
        self.radius : np.array
        self.cog : np.array
        self.minCoords : np.array
        self.maxCoords : np.array
        self.resolution : np.array
        self.deltas : np.array
        self.set_box_properties()

        if VERBOSE: print(
            f"...... Resolution = {volpot.format_vector_str(self.resolution)}, Deltas = {volpot.format_vector_str(self.deltas)}, "\
            f"COG = {volpot.format_vector_str(self.cog)}, radius = {self.radius:.3f}", flush = True
        )


        ###############################
        grid_size = np.prod(self.resolution)
        if grid_size > WARNING_GRID_SIZE:
            print()
            while True:
                choice = input(f">>> WARNING: resulting grid would contain {grid_size/1e6:.2f} million points. Proceed? [Y/N]\n").upper()
                if choice == 'Y': break
                if choice == 'N': exit()

        ###############################
        self.relevant_atoms = mda.core.groups.AtomGroup
        self.relevant_atoms_broad = mda.core.groups.AtomGroup
        self.select_relevant_atoms()

        self.trimming_mask : np.array
        self.coords = self.minCoords + np.indices(self.resolution).T * self.deltas
        self.create_trimming_mask(self.trimming_dist)


    def create_trimming_mask(self, min_dist):
        timer = volpot.Timer("...>>> MS: Generating mask...")
        self.trimming_mask = np.zeros(self.resolution, dtype = bool)
        if DO_TRIMMING_SPHERE:    self.trim_sphere()
        if DO_TRIMMING_OCCUPANCY: self.trim_occupancy(min_dist)
        if DO_TRIMMING_RNDS:      self.trim_rnds()
        timer.end()

        if not SAVE_CACHED_MASK: return

        if VERBOSE: timer = volpot.Timer("......   Saving mask cache...")
        if DO_DX_OUTPUT:
            path_cached_mask = self.folder_potentials / f"{self.name}.mask.dx"
            self.write_dx(path_cached_mask, self.grid)
        if DO_MRC_OUTPUT:
            path_cached_mask = self.folder_potentials / f"{self.name}.mask.mrc"
            self.write_mrc(path_cached_mask, self.grid)
        if VERBOSE: timer.end()


    def trim_occupancy(self, radius):
        sk = volpot.SphereKernel(radius, self.deltas, bool)
        sk.link_to_grid(self.trimming_mask, self.minCoords)
        for a in self.relevant_atoms_broad:
            sk.stamp(a.position)


    ######################### SPECIFIC METHODS (OVERRIDE TO USE)
    def set_box_properties(self): return
    def select_relevant_atoms(self): return
    def trim_sphere(self): return
    def trim_rnds(self): return


    def write_dx(self, path_dx, grid_data):
        ########### init values
        xres,yres,zres = self.resolution
        xmin,ymin,zmin = self.minCoords
        dx,dy,dz = self.deltas

        if grid_data.dtype == np.float32:
            dtype = '"float"'
            fmt = "%.3f"
        elif grid_data.dtype == int:
            dtype = '"int"'
            fmt = "%i"

        header = self.dx_comments + '\n'.join((
            f"object 1 class gridpositions counts  {xres} {yres} {zres}",
            f"origin {xmin:6f} {ymin:6f} {zmin:6f}",
            f"delta  {dx} 0 0",
            f"delta  0 {dy} 0",
            f"delta  0 0 {dz}",
            f"object 2 class gridconnections counts  {xres} {yres} {zres}",
            f"object 3 class array type {dtype} rank 0 items {np.prod(self.resolution)} data follows",
        ))

        ########### reshape the grid array
        grid_size = np.prod(grid_data.shape)
        dx_rows = grid_size // 3

        truncated_arr, extra_arr = np.split(grid_data.flatten(), [3*dx_rows])
        data_out = truncated_arr.reshape(dx_rows, 3)
        last_row = extra_arr.reshape(1, len(extra_arr))

        ########### export reshaped data
        with open(path_dx, "wb") as file:
            np.savetxt(
                file, data_out, fmt = fmt, delimiter = '\t',
                header = header, comments = ''
            )
            np.savetxt(
                file, last_row, fmt = fmt, delimiter = '\t',
                footer = self.dx_footer, comments = ''
            )


    # --------------------------------------------------------------------------
    def write_mrc(self, path_mrc, grid_data):
        with mrc.mrcfile.new(path_mrc, overwrite = True) as mrc_file:
            mrc_file.set_data(grid_data.T.astype(np.float32))
            mrc_file.voxel_size = self.deltas.tolist()
            mrc_file.header["origin"]['x'] = self.minCoords[0]
            mrc_file.header["origin"]['y'] = self.minCoords[1]
            mrc_file.header["origin"]['z'] = self.minCoords[2]
            mrc_file.update_header_from_data()
            mrc_file.update_header_stats()



# //////////////////////////////////////////////////////////////////////////////
class MS_PocketSphere(MolecularSystem):
    def set_box_properties(self):
        self.cog = np.array([self.metadata["xcog"], self.metadata["ycog"], self.metadata["zcog"]])
        self.radius = self.metadata["radius"]
        self.minCoords = self.cog - self.radius
        self.maxCoords = self.cog + self.radius
        box_size = self.maxCoords - self.minCoords

        if self.metadata["default_res"]:
            self.resolution = np.array([GRID_XRES_PS, GRID_YRES_PS, GRID_ZRES_PS])
            self.deltas = box_size / self.resolution
        else:
            self.deltas = np.array([GRID_DX_PS, GRID_DY_PS, GRID_DZ_PS])
            self.resolution = np.round(box_size / self.deltas).astype(int)

        print("...>>> MS: PocketSphere mode. Info:", flush = True)


    def select_relevant_atoms(self):
        xcog, ycog, zcog = self.cog
        self.relevant_atoms = self.system.select_atoms(
            f"{self.macro_query} and point {xcog} {ycog} {zcog} {self.radius}"
        )
        self.relevant_atoms_broad = self.system.select_atoms(
            f"{self.macro_query} and point {xcog} {ycog} {zcog} {self.radius + self.trimming_dist}"
        )

    def trim_sphere(self):
        shifted_coords = self.coords - self.cog
        dist_from_cog = volpot.get_norm(shifted_coords)
        self.trimming_mask[dist_from_cog > self.radius] = True

    def trim_rnds(self):
        """perform a random search to remove isolated regions"""
        visited = np.zeros(self.resolution, dtype = bool)

        directions = np.array([[i,j,k] for i in range(-1,2) for j in range(-1,2) for k in range(-1,2) if i&j&k])

        xres, yres, zres = self.resolution
        xcog, ycog, zcog = np.floor(self.resolution / 2).astype(int)
        cog_cube = set((x,y,z)
            for x in range(xcog - COG_CUBE_RADIUS, xcog + COG_CUBE_RADIUS + 1)
            for y in range(ycog - COG_CUBE_RADIUS, ycog + COG_CUBE_RADIUS + 1)
            for z in range(zcog - COG_CUBE_RADIUS, zcog + COG_CUBE_RADIUS + 1)
        )
        queue = cog_cube.copy()

        search_dist = np.full(self.resolution, np.inf)
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
                if search_dist[neigh] > MAX_RNDS_DIST: continue
                if visited[neigh]: continue
                if self.trimming_mask[neigh]: continue

                queue.add(neigh)

        self.trimming_mask[np.logical_not(visited)] = True


# //////////////////////////////////////////////////////////////////////////////
class MS_Whole(MolecularSystem):
    def set_box_properties(self):
        self.maxCoords = np.max(self.system.coord.positions, axis = 0) + EXTRA_BOX_SIZE
        self.minCoords = np.min(self.system.coord.positions, axis = 0) - EXTRA_BOX_SIZE
        self.radius = volpot.get_norm(self.maxCoords - self.minCoords) / 2
        self.cog = (self.maxCoords + self.minCoords) / 2
        box_size = self.maxCoords - self.minCoords

        if self.metadata["default_res"]:
            self.resolution = np.array([GRID_XRES_WM, GRID_YRES_WM, GRID_ZRES_WM])
            self.deltas = box_size / self.resolution
        else:
            self.deltas = np.array([GRID_DX_WM, GRID_DY_WM, GRID_DZ_WM])
            self.resolution = np.round(box_size / self.deltas).astype(int)

        print("...>>> MS: Whole mode. Info:", flush = True)

    def select_relevant_atoms(self):
        self.relevant_atoms = self.system.select_atoms(self.macro_query)
        self.relevant_atoms_broad = self.relevant_atoms


################################################################################
