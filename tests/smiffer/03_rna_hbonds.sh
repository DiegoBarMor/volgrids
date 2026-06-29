#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 3: Benchmark RNA, Hydrogen bonds details"

fapbs="testdata/smiffer/apbs"
fpdb_orig="testdata/smiffer/pdb_orig"
fpdb_clean="testdata/smiffer/pdb_clean"


fout="testdata/smiffer/whole-hbond_details"
mkdir -p $fout

conf_just_hbond="OUT_FORMAT=CMAP SMIF_STK=false SMIF_HBA=true SMIF_HBD=true SMIF_HPHOB=false SMIF_HPHIL=false SMIF_APBS=false"

names=(1akx 1i9v 2esj 4f8u 5bjo 5kx9 6tf3 7oax0 8eyv)
for name in "${names[@]}"; do
    cp "$fpdb_orig/$name.pdb" "$fout/$name.pdb"

    ##### PART 0: H-bonds for only nucleobases, all residues
    python3 volgrids smiffer "$fpdb_clean/$name.pdb" -o "$fout" \
        --pack -c "$conf_just_hbond" SMIF_HB_ONLY_NBASE=true
    mv "$fout/$name.all.cmap" "$fout/$name.nbases.cmap"


    ##### PART 1: H-bonds for only nucleobases, non-base-paired residues only
    if ! residues="$(python3 volgrids smutils res_nobp "$fpdb_clean/$name.pdb" 2>&1)"; then
        rc=$?
        printf 'res_nobp failed (exit %s). Output:\n%s\n' "$rc" "$residues" >&2
        exit $rc
    fi
    echo "... non-base-paired residues for $name: $residues"

    python3 volgrids smiffer "$fpdb_clean/$name.pdb" -o "$fout" \
        --pack -c "$conf_just_hbond" SMIF_HB_ONLY_NBASE=true -r "$residues"
    mv "$fout/$name.all.cmap" "$fout/$name.nbases.nobp.cmap"


    ##### PART 2: H-bonds for all residues + APBS
    python3 volgrids smiffer "$fpdb_clean/$name.pdb" -o "$fout" \
        -c "$conf_just_hbond" SMIF_APBS=true -a "$fapbs/$name.pdb.mrc"
done
