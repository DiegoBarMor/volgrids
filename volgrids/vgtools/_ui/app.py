import volgrids as vg
import volgrids.vgtools as vgt

# //////////////////////////////////////////////////////////////////////////////
class AppVGTools(vg.App):
    CONFIG_MODULES = (vg,)
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

        if vgt.OPERATION == "average":
            print(f">>> Averaging CMAP file: {vgt.PATH_AVERAGE_IN}")
            vgt.VGOperations.average(vgt.PATH_AVERAGE_IN, vgt.PATH_AVERAGE_OUT)
            return

        if vgt.OPERATION == "summary":
            print(f">>> Grid summary: {vgt.PATH_SUMMARY_IN}")
            vgt.VGOperations.summary(vgt.PATH_SUMMARY_IN)
            return

        if vgt.OPERATION == "compare":
            print(f">>> Comparing grids: {vgt.PATH_COMPARE_IN_0} vs {vgt.PATH_COMPARE_IN_1} (threshold={vgt.THRESHOLD_COMPARE:2.2e})")
            result = vgt.VGOperations.compare(vgt.PATH_COMPARE_IN_0, vgt.PATH_COMPARE_IN_1, vgt.THRESHOLD_COMPARE)

            for message in result.messages:
                print(f"...>>> {message}")
            if result.npoints_total == 0: return

            print(
                f"...>>> {result.npoints_diff}/{result.npoints_total} points differ " +\
                f"({100 * result.npoints_diff / result.npoints_total:.2f}%)\n" +\
                f"...>>> Accumulated difference: {result.cumulative_diff:2.2e} " +\
                f"(avg {result.avg_diff:2.2e} per point)"
            )
            return

        if vgt.OPERATION == "rotate":
            print(f">>> Rotating grid: {vgt.PATH_ROTATE_IN} by {vgt.ROTATE_XY}° (xy), {vgt.ROTATE_YZ}° (yz), {vgt.ROTATE_XZ}° (xz)")
            vgt.VGOperations.rotate(vgt.PATH_ROTATE_IN, vgt.PATH_ROTATE_OUT, vgt.ROTATE_XY, vgt.ROTATE_YZ, vgt.ROTATE_XZ)
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
