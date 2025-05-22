import volpot as vp

# ------------------------------------------------------------------------------
def open_grid(path) -> tuple[str, vp.Grid]:
    if path.suffix == ".dx":
        return "DX", vp.read_dx(path)

    if path.suffix == ".mrc":
        return "MRC", vp.read_mrc(path)

    elif path.suffix == ".cmap":
        keys = vp.get_cmap_keys(path)
        if not keys: raise ValueError(f"Empty cmap file: {path}")
        return "CMAP", vp.read_cmap(path, keys[0])

    raise ValueError(f"Unrecognized file format: {path.suffix}")


################################################################################
if __name__ == "__main__":
    args = vp.args_tools()

    # --------------------------------------------------------------------------
    if str(args.to_mrc) != '.':
        path_in = args.to_mrc
        path_out = path_in.with_suffix(".mrc")
        kind, grid = open_grid(path_in)
        print(f">>> Converting {kind} file into MRC: {path_in}")
        vp.write_mrc(path_out, grid)


    # --------------------------------------------------------------------------
    if str(args.to_cmap) != '.':
        path_in = args.to_cmap
        path_out = path_in.with_suffix(".cmap")
        print(">>> Converting DX file to CMAP:", path_in)
        vp.write_cmap(path_out, grid, path_in.stem)


    # --------------------------------------------------------------------------
    if isinstance(args.pack, list):
        if len(args.pack) < 2:
            raise ValueError("Packing requires at least two files: one for the output and one or more for the input.")

        it = iter(args.pack)
        path_out = next(it)

        print(f">>> Packing {[p.name for p in args.pack[1:]]} into '{path_out}'")
        for path_in in it:
            if not path_in.exists():
                print("...>>> Skipping non-existing file:", path_in)
                continue
            kind, grid = open_grid(path_in)
            vp.write_cmap(path_out, grid, path_in.stem)


    # --------------------------------------------------------------------------
    if str(args.unpack) != '.':
        path_in = args.unpack

        keys = vp.get_cmap_keys(path_in)
        print(f">>> Unpacking {keys} from '{path_in}' into '{path_out}'")
        for key in keys:
            path_out = path_in.parent / f"{key}.cmap"
            grid = vp.read_cmap(path_in, key)
            vp.write_cmap(path_out, grid)


################################################################################
