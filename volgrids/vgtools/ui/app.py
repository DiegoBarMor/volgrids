import volgrids.vgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.App):
    _CLASS_PARAM_HANDLER = vgt.ParamHandlerVGTools

    # --------------------------------------------------------------------------
    def run(self) -> None:
        if vgt.OPERATION == "convert":
            self._run_convert()
            return

        if vgt.OPERATION == "pack":
            print(f">>> Packing {len(vgt.PATHS_PACK_IN)} grids into '{vgt.PATH_PACK_OUT}'")
            vgt.VGOperations.pack(vgt.PATHS_PACK_IN, vgt.PATH_PACK_OUT)
            return

        if vgt.OPERATION == "unpack":
            print(f">>> Unpacking '{vgt.PATH_UNPACK_IN}' into '{vgt.PATH_UNPACK_OUT}'")
            vgt.VGOperations.unpack(vgt.PATH_UNPACK_IN, vgt.PATH_UNPACK_OUT)
            return

        if vgt.OPERATION == "fix_cmap":
            print(f">>> Fixing CMAP file: {vgt.PATH_FIXCMAP_IN}")
            vgt.VGOperations.fix_cmap(vgt.PATH_FIXCMAP_IN, vgt.PATH_FIXCMAP_OUT)
            return

        raise ValueError(f"Unknown mode: {vgt.OPERATION}")


    # --------------------------------------------------------------------------
    def _run_convert(self):
        def _convert(path_out, fmt_out: vg.GridFormat):
            if path_out is None: return
            print(f">>> Converting {vgt.PATH_CONVERT_IN} file to {fmt_out.name}: {path_out}")
            vgt.VGOperations.convert(vgt.PATH_CONVERT_IN, path_out, fmt_out)

        _convert(vgt.PATH_CONVERT_DX,   vg.GridFormat.DX)
        _convert(vgt.PATH_CONVERT_MRC,  vg.GridFormat.MRC)
        _convert(vgt.PATH_CONVERT_CCP4, vg.GridFormat.CCP4)
        _convert(vgt.PATH_CONVERT_CMAP, vg.GridFormat.CMAP)


# //////////////////////////////////////////////////////////////////////////////
