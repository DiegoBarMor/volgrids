#!/bin/bash
set -euo pipefail

config="examples/interfaces/config.ini"
fdata="testdata/smiffer/interfaces/prot_rna"
fp2n=$fdata/prot2nucl
fn2p=$fdata/nucl2prot
mkdir -p $fp2n $fn2p


######### PART 1: calculate SMIF and OG grids
python3 volgrids smiffer $fdata/prot.pdb -c "$config" -o $fp2n
python3 volgrids smiffer $fdata/rna.pdb  -c "$config" -o $fn2p

python3 volgrids smutils occupancy $fdata/prot.pdb -c "$config" -o $fn2p
python3 volgrids smutils occupancy $fdata/rna.pdb  -c "$config" -o $fp2n


######### PART 2: overlap the SMIFs with the OGs
overlap() {
    python3 volgrids vgtools op mul "$1" "$2" "$3"
}
overlap $fp2n/prot.stk.smif.bin $fp2n/rna.stk.og.bin $fp2n/p2n.stk.bin
overlap $fp2n/prot.hba.smif.bin $fp2n/rna.hba.og.bin $fp2n/p2n.hbd.bin
overlap $fp2n/prot.hbd.smif.bin $fp2n/rna.hbd.og.bin $fp2n/p2n.hba.bin

overlap $fn2p/rna.stk.smif.bin $fn2p/prot.stk.og.bin $fn2p/n2p.stk.bin
overlap $fn2p/rna.hba.smif.bin $fn2p/prot.hba.og.bin $fn2p/n2p.hbd.bin
overlap $fn2p/rna.hbd.smif.bin $fn2p/prot.hbd.og.bin $fn2p/n2p.hba.bin


######### PART 3: Pack the BIN files into CMAP
pack() {
    python3 volgrids vgtools pack "$1" "$2" "$3" -o "$4"
}
pack $fp2n/prot.stk.smif.bin $fp2n/prot.hbd.smif.bin $fp2n/prot.hba.smif.bin $fdata/prot.smif.cmap
pack $fp2n/rna.stk.og.bin    $fp2n/rna.hba.og.bin    $fp2n/rna.hbd.og.bin    $fdata/rna.og.cmap
pack $fn2p/rna.stk.smif.bin  $fn2p/rna.hbd.smif.bin  $fn2p/rna.hba.smif.bin  $fdata/rna.smif.cmap
pack $fn2p/prot.stk.og.bin   $fn2p/prot.hba.og.bin   $fn2p/prot.hbd.og.bin   $fdata/prot.og.cmap
pack $fp2n/p2n.stk.bin $fp2n/p2n.hbd.bin $fp2n/p2n.hba.bin $fdata/prot_smif.rna_og.cmap
pack $fn2p/n2p.stk.bin $fn2p/n2p.hbd.bin $fn2p/n2p.hba.bin $fdata/rna_smif.prot_og.cmap


######### PART 4: Calculate segmentation (from "examples/bin_format") for one of the grids
path_smif="$fp2n/prot.stk.smif.bin"
path_clusters="$fdata/prot.stk.clusters.bin"
isovalue=0.1

python3 volgrids vgtools segment "$path_smif" "$path_clusters" -i $isovalue
python3 examples/bin_format/expand_mrc2cmap.py "$path_clusters" ### [WIP] move this to a proper operation


######### CLEANUP
rm -rf $fp2n $fn2p
