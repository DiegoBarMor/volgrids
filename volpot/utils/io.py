import json

################################################################################
def read_json(path_json):
    with open(path_json, 'r') as file:
        return json.load(file)

# ------------------------------------------------------------------------------
def write_json(path_json, content):
    with open(path_json, 'w') as file:
        file.write(json.dumps(content))

# ------------------------------------------------------------------------------
def save_metadata(metadata):
    meta = metadata.copy()
    for k in ("pdb", "out", "apbs", "meta"):
        meta[k] = str(meta[k])
    write_json(meta["meta"], meta)


################################################################################
