import volgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class VGTools:
    def __init__(self):
        self.meta = vgt.VGToolsArgsParser()


    # --------------------------------------------------------------------------
    def run(self):
        if self.meta.mode == "convert":
            self._run_convert()
            return

        if self.meta.mode == "pack":
            self._run_pack()
            return

        if self.meta.mode == "unpack":
            self._run_unpack()
            return

        if self.meta.mode == "fix-cmap":
            self._run_fix_cmap()
            return

        raise ValueError(f"Unknown mode: {self.meta.mode}")


    # --------------------------------------------------------------------------
    def _run_convert(self):
        kind, grid = vg.GridIO.read_auto(self.meta.path_in)

        if self.meta.path_dx is not None:
            print(f">>> Converting {kind.name} file to DX: {self.meta.path_dx}")
            vg.GridIO.write_dx(self.meta.path_dx, grid)

        if self.meta.path_mrc is not None:
            print(f">>> Converting {kind.name} file to MRC: {self.meta.path_mrc}")
            vg.GridIO.write_mrc(self.meta.path_mrc, grid)

        if self.meta.path_cmap is not None:
            print(f">>> Converting {kind.name} file to CMAP: {self.meta.path_cmap}")
            vg.GridIO.write_cmap(self.meta.path_cmap, grid, self.meta.path_in.stem)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        paths_in = self.meta.paths_pack_in
        path_out = self.meta.path_pack_out

        resolution = None
        print(f">>> Packing {len(paths_in)} grids into '{path_out}'")
        for path_in in paths_in:
            if not path_in.exists():
                print("...>>> Skipping non-existing file:", path_in)
                continue

            _, grid = vg.GridIO.read_auto(path_in)
            if resolution is None:
                resolution = (grid.xres, grid.yres, grid.zres)

            new_res = (grid.xres, grid.yres, grid.zres)
            if new_res != resolution:
                print(
                    f">>> Warning: Grid {path_in} has different resolution {new_res} than the first grid {resolution}. " +\
                    "Chimera won't recognize it as a volume series and open every grid in a separate representation." +\
                    "Use `smiffer.py fix-cmap` if you want to fix this."
                )
            key = str(path_in.parent / path_in.stem).replace(' ', '_').replace('/', '_').replace('\\', '_')
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        path_in = self.meta.path_unpack_in
        folder_out = self.meta.path_unpack_out

        keys = vg.GridIO.get_cmap_keys(path_in)
        print(f">>> Unpacking {keys} from '{path_in}' into '{folder_out}'")
        for key in keys:
            path_out = folder_out / f"{key}.cmap"
            grid = vg.GridIO.read_cmap(path_in, key)
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    def _run_fix_cmap(self):
        path_in = self.meta.path_fixcmap_in
        path_out = self.meta.path_fixcmap_out

        resolution = None
        keys = vg.GridIO.get_cmap_keys(path_in)
        print(f">>> Fixing CMAP file: {path_in}")
        for key in keys:
            grid = vg.GridIO.read_cmap(path_in, key)

            minCoords = (grid.xmin, grid.ymin, grid.zmin)
            maxCoords = (grid.xmax, grid.ymax, grid.zmax)
            if resolution is None:
                resolution = (grid.xres, grid.yres, grid.zres)

            grid.reshape(minCoords, maxCoords, resolution)
            vg.GridIO.write_cmap(path_out, grid, key)


# //////////////////////////////////////////////////////////////////////////////
