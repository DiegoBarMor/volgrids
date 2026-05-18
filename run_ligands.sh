#!/bin/bash
set -euo pipefail

data=...
rm -f $data/*.chem $data/*.mrc

cd $data
for name in $(ls ./*.pdb); do
    # clear
    echo $name
    stem=$(basename $name .pdb)

    volgrids smutils chemgen $name 2>/dev/null
    if [[ "$name" == *_ARG_* ]]; then
        echo Manually adding flat group for arginine...
        echo "ARG: CZ NE NH1 NH2" >> $stem.chem
    fi

    volgrids smiffer ligand $name -b $stem.chem -c config.ini

#     cat >tmp.vmd <<EOF
# mol new $stem.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
# mol addfile $stem.stacking.mrc type ccp4 first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
# mol delrep 0 top
# mol representation Licorice 0.300000 12.000000 12.000000
# mol color Name
# mol selection {all}
# mol material Opaque
# mol addrep top
# mol representation Isosurface 1.000000 0 2 0 1 1
# mol color Name
# mol selection {all}
# mol material Transparent
# mol addrep top
# EOF
#     vmd -e tmp.vmd >/dev/null
done

# rm -f tmp.vmd
