import volgrids as vg
import volgrids.smutils as su

# //////////////////////////////////////////////////////////////////////////////
class AppSMUtils(vg.App):
    CONFIG_MODULES = (vg,)
    _CLASS_PARAM_HANDLER = su.ParamHandlerSMUtils

    # --------------------------------------------------------------------------
    def run(self) -> None:
        if su.OPERATION == "resids_nonbp":
            su.SMOperations.print_resids_nonbp(su.PATH_STRUCT)
            return

        raise ValueError(f"Unknown mode: {su.OPERATION}")


# //////////////////////////////////////////////////////////////////////////////
