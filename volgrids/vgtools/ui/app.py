import volgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.App):
    _CLASS_PARAM_HANDLER = vgt.ParamHandlerVGTools

    # --------------------------------------------------------------------------
    def run(self):
        mode = vg.USER_MODE.lower()
        if mode == "convert":
            self._run_convert()
            return

        if mode == "pack":
            self._run_pack()
            return

        if mode == "unpack":
            self._run_unpack()
            return

        if mode == "fix_cmap":
            self._run_fix_cmap()
            return

        raise ValueError(f"Unknown mode: {mode}")


    # --------------------------------------------------------------------------
    def _run_convert(self):
        kind, grid = vg.GridIO.read_auto(vgt.PATH_CONVERT_IN)

        if vgt.PATH_CONVERT_DX is not None:
            print(f">>> Converting {kind.name} file to DX: {vgt.PATH_CONVERT_DX}")
            vg.GridIO.write_dx(vgt.PATH_CONVERT_DX, grid)

        if vgt.PATH_CONVERT_MRC is not None:
            print(f">>> Converting {kind.name} file to MRC: {vgt.PATH_CONVERT_MRC}")
            vg.GridIO.write_mrc(vgt.PATH_CONVERT_MRC, grid)

        if vgt.PATH_CONVERT_CCP4 is not None:
            print(f">>> Converting {kind.name} file to CCP4: {vgt.PATH_CONVERT_CCP4}")
            vg.GridIO.write_ccp4(vgt.PATH_CONVERT_CCP4, grid)

        if vgt.PATH_CONVERT_CMAP is not None:
            print(f">>> Converting {kind.name} file to CMAP: {vgt.PATH_CONVERT_CMAP}")
            vg.GridIO.write_cmap(vgt.PATH_CONVERT_CMAP, grid, vgt.PATH_CONVERT_IN.stem)


    # --------------------------------------------------------------------------
    def _run_pack(self):
        paths_in = vgt.PATHS_PACK_IN
        path_out = vgt.PATH_PACK_OUT

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
                    "Use `smiffer.py fix_cmap` if you want to fix this."
                )
            key = str(path_in.parent / path_in.stem).replace(' ', '_').replace('/', '_').replace('\\', '_')
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    def _run_unpack(self):
        path_in    = vgt.PATH_UNPACK_IN
        folder_out = vgt.PATH_UNPACK_OUT

        keys = vg.GridIO.get_cmap_keys(path_in)
        print(f">>> Unpacking {keys} from '{path_in}' into '{folder_out}'")
        for key in keys:
            path_out = folder_out / f"{key}.cmap"
            grid = vg.GridIO.read_cmap(path_in, key)
            vg.GridIO.write_cmap(path_out, grid, key)


    # --------------------------------------------------------------------------
    def _run_fix_cmap(self):
        path_in  = vgt.PATH_FIXCMAP_IN
        path_out = vgt.PATH_FIXCMAP_OUT

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
