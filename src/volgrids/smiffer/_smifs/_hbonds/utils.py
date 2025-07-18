import volgrids as vg

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
def safe_return_coords(atomgroup, sel_string):
    atoms = atomgroup.select_atoms(sel_string)
    if len(atoms) == 0: return None
    return atoms.center_of_geometry()


# ------------------------------------------------------------------------------
def get_direction_vector(vec_origin, vec_head):
    return vg.Math.normalize(vec_head - vec_origin)


# ------------------------------------------------------------------------------
