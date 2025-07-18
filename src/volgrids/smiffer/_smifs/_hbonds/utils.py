# ------------------------------------------------------------------------------
def str_prev_residue(res):
    return f"segid {res.segid} and resid {res.resid - 1}"


# ------------------------------------------------------------------------------
def str_this_residue(res):
    return f"segid {res.segid} and resid {res.resid} and resname {res.resname}"


# ------------------------------------------------------------------------------
def str_next_residue(res):
    return f"segid {res.segid} and resid {res.resid + 1}"


# ------------------------------------------------------------------------------
