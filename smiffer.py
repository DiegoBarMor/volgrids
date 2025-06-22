import volgrids as vg
import volgrids.smiffer as sm

################################################################################
if __name__ == "__main__":
    meta = sm.SmifferArgsParser()

    for name, value in meta.debug_vars.items():
        if hasattr(vg, name):
            setattr(vg, name, value)
        elif hasattr(sm, name):
            setattr(sm, name, value)
        else:
            print(f"Warning: {name} is not a valid volgrids variable. Skipping.")
            continue

    sm.SmifferCalculator(meta).run()


################################################################################
