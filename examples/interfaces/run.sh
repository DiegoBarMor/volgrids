#!/bin/bash
set -euo pipefail

config="examples/interfaces/config.ini"
fdata="testdata/smiffer/interfaces/prot_rna"
fp2n=$fdata/prot2nucl
fn2p=$fdata/nucl2prot

mkdir -p $fp2n $fn2p

python3 volgrids smiffer $fdata/prot.pdb -c "$config" -o $fp2n
python3 volgrids smiffer $fdata/rna.pdb  -c "$config" -o $fn2p

python3 volgrids smutils occupancy $fdata/prot.pdb -c "$config" -o $fn2p
python3 volgrids smutils occupancy $fdata/rna.pdb  -c "$config" -o $fp2n

overlap() {
    python3 volgrids vgtools op mul "$1" "$2" "$3"
}
pack() {
    python3 volgrids vgtools pack "$1" "$2" "$3" -o "$4"
}

overlap $fp2n/prot.stacking.mrc    $fp2n/rna.og_stacking.mrc    $fp2n/p2n.stk.mrc
overlap $fp2n/prot.hbdonors.mrc    $fp2n/rna.og_hbacceptors.mrc $fp2n/p2n.hbd.mrc
overlap $fp2n/prot.hbacceptors.mrc $fp2n/rna.og_hbdonors.mrc    $fp2n/p2n.hba.mrc

overlap $fn2p/rna.stacking.mrc    $fn2p/prot.og_stacking.mrc    $fn2p/n2p.stk.mrc
overlap $fn2p/rna.hbdonors.mrc    $fn2p/prot.og_hbacceptors.mrc $fn2p/n2p.hbd.mrc
overlap $fn2p/rna.hbacceptors.mrc $fn2p/prot.og_hbdonors.mrc    $fn2p/n2p.hba.mrc

pack $fp2n/prot.stacking.mrc $fp2n/prot.hbdonors.mrc $fp2n/prot.hbacceptors.mrc $fdata/prot.smif.cmap
pack $fp2n/rna.og_stacking.mrc $fp2n/rna.og_hbacceptors.mrc $fp2n/rna.og_hbdonors.mrc $fdata/rna.og.cmap
pack $fn2p/rna.stacking.mrc $fn2p/rna.hbdonors.mrc $fn2p/rna.hbacceptors.mrc $fdata/rna.smif.cmap
pack $fn2p/prot.og_stacking.mrc $fn2p/prot.og_hbacceptors.mrc $fn2p/prot.og_hbdonors.mrc $fdata/prot.og.cmap
pack $fp2n/p2n.stk.mrc $fp2n/p2n.hbd.mrc $fp2n/p2n.hba.mrc $fdata/prot_smif.rna_og.cmap
pack $fn2p/n2p.stk.mrc $fn2p/n2p.hbd.mrc $fn2p/n2p.hba.mrc $fdata/rna_smif.prot_og.cmap

rm -rf $fp2n $fn2p
