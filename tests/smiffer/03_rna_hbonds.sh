#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 3: Benchmark RNA, Hydrogen bonds details"

fapbs="testdata/smiffer/apbs"
fpdb_orig="testdata/smiffer/pdb_orig"
fpdb_clean="testdata/smiffer/pdb_clean"


fout="testdata/smiffer/whole-hbond_details"
mkdir -p $fout

conf_just_hbond="DO_SMIF_STACKING=false DO_SMIF_HBA=true DO_SMIF_HBD=true DO_SMIF_HYDROPHOBIC=false DO_SMIF_HYDROPHILIC=false DO_SMIF_APBS=false"

names=(1akx 1i9v 2esj 4f8u 5bjo 5kx9 6tf3 7oax0 8eyv)
for name in "${names[@]}"; do
    cp "$fpdb_orig/$name.pdb" "$fout/$name.pdb"

    ##### PART 0: H-bonds for only nucleobases, all residues
    python3 volgrids smiffer rna  "$fpdb_clean/$name.pdb" -o "$fout" \
        -c "$conf_just_hbond" HBONDS_ONLY_NUCLEOBASE=true
    mv "$fout/$name.cmap" "$fout/$name.nbases.cmap"


    ##### PART 1: H-bonds for only nucleobases, non-base-paired residues only
    if ! idxs="$(python3 volgrids smutils resids_nonbp "$fpdb_clean/$name.pdb" 2>&1)"; then
        rc=$?
        printf 'resids_nonbp failed (exit %s). Output:\n%s\n' "$rc" "$idxs" >&2
        exit $rc
    fi
    echo "... non-base-paired residue indices for $name: $idxs"

    python3 volgrids smiffer rna  "$fpdb_clean/$name.pdb" -o "$fout" \
        -c "$conf_just_hbond" HBONDS_ONLY_NUCLEOBASE=true -i "$idxs"
    mv "$fout/$name.cmap" "$fout/$name.nbases.nonbp.cmap"


    ##### PART 2: H-bonds for all residues + APBS
    python3 volgrids smiffer rna  "$fpdb_clean/$name.pdb" -o "$fout" \
        -c "$conf_just_hbond" DO_SMIF_APBS=true -a "$fapbs/$name.pdb.mrc"
done
