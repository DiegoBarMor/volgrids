#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 5: Cavities"

fpdb="testdata/smiffer/pdb-nosolv"
fpdb_orig="testdata/_raw_input/smiffer_benchmark"
fout="testdata/smiffer/cavities"
fout_benchmark=$fout/benchmark
fout_pocketsphere=$fout/pocket_sphere
fout_options=$fout/options
rm -rf $fout
mkdir -p $fout_benchmark $fout_pocketsphere $fout_options


############################# BENCHMARK SYSTEMS
tmp_config_benchmark="$fout_benchmark/config.tmp"

cat > $tmp_config_benchmark <<- EOM
DO_SMIF_STACKING=False
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False

DO_TRIMMING_OCCUPANCY=True
DO_TRIMMING_CAVITIES=False
SAVE_TRIMMING_MASK=True
SAVE_CAVITIES=True
TRIMMING_CAVITIES_THRESHOLD=3
EOM

for name in 1bg0 1eby 1ehe 1h7l 1iqj 1ofz 3dd0 3ee4 5m9w 6e9a ; do
    python3 volgrids smiffer prot  $fpdb/$name.pdb -o $fout_benchmark --config $tmp_config_benchmark
    cp $fpdb_orig/$name.pdb $fout_benchmark/
done
for name in 1akx 1i9v 2esj 4f8u 5bjo 5kx9 6tf3 7oax0 7oax1 8eyv; do
    python3 volgrids smiffer rna  $fpdb/$name.pdb -o $fout_benchmark --config $tmp_config_benchmark
    cp $fpdb_orig/$name.pdb $fout_benchmark/
done

rm -f $tmp_config_benchmark


############################# OPTIONS
write_config() {
    local cfg_path="$1"
    local do_trimming_occupancy="$2"
    local do_trimming_cavities="$3"
    local save_trimming_mask="$4"
    local save_cavities="$5"
    local trimming_threshold="$6"

    cat > "$cfg_path" <<- EOM
DO_SMIF_STACKING=True
DO_SMIF_HBA=False
DO_SMIF_HBD=False
DO_SMIF_HYDROPHOBIC=False
DO_SMIF_HYDROPHILIC=False
DO_SMIF_APBS=False

### configured trimming/cavity options
DO_TRIMMING_OCCUPANCY=${do_trimming_occupancy}
DO_TRIMMING_CAVITIES=${do_trimming_cavities}
SAVE_TRIMMING_MASK=${save_trimming_mask}
SAVE_CAVITIES=${save_cavities}
TRIMMING_CAVITIES_THRESHOLD=${trimming_threshold}
EOM
}

path_pdb="$fpdb/1iqj.pdb"


######################### booleans
tmp_config_options_00="$fout_options/config_00.tmp"
tmp_config_options_01="$fout_options/config_01.tmp"
tmp_config_options_02="$fout_options/config_02.tmp"
tmp_config_options_03="$fout_options/config_03.tmp"
tmp_config_options_04="$fout_options/config_04.tmp"
write_config "$tmp_config_options_00" False True  True True  3 # DO_TRIMMING_OCCUPANCY=false, so cavities should be ignored
write_config "$tmp_config_options_01" True  False True False 3 # cavities shouldn't be saved and shouldn't affect the trimming
write_config "$tmp_config_options_02" True  True  True False 3 # cavities shouldn't be saved but should affect the trimming
write_config "$tmp_config_options_03" True  False True True  3 # cavities should be saved but shouldn't affect the trimming
write_config "$tmp_config_options_04" True  True  True True  3 # cavities should be saved and should affect the trimming

cp $path_pdb $fout_options/config_00.pdb
cp $path_pdb $fout_options/config_01.pdb
cp $path_pdb $fout_options/config_02.pdb
cp $path_pdb $fout_options/config_03.pdb
cp $path_pdb $fout_options/config_04.pdb
python3 volgrids smiffer prot $fout_options/config_00.pdb --config $tmp_config_options_00
python3 volgrids smiffer prot $fout_options/config_01.pdb --config $tmp_config_options_01
python3 volgrids smiffer prot $fout_options/config_02.pdb --config $tmp_config_options_02
python3 volgrids smiffer prot $fout_options/config_03.pdb --config $tmp_config_options_03
python3 volgrids smiffer prot $fout_options/config_04.pdb --config $tmp_config_options_04

rm -f $tmp_config_options_00 $tmp_config_options_01 $tmp_config_options_02 $tmp_config_options_03 $tmp_config_options_04


######################### threshold
tmp_config_options_t0="$fout_options/config_t0.tmp"
tmp_config_options_t1="$fout_options/config_t1.tmp"
tmp_config_options_t2="$fout_options/config_t2.tmp"
tmp_config_options_t3="$fout_options/config_t3.tmp"
tmp_config_options_t4="$fout_options/config_t4.tmp"
write_config "$tmp_config_options_t0" True True True True 0
write_config "$tmp_config_options_t1" True True True True 1
write_config "$tmp_config_options_t2" True True True True 2
write_config "$tmp_config_options_t3" True True True True 3
write_config "$tmp_config_options_t4" True True True True 4

cp $path_pdb $fout_options/config_t0.pdb
cp $path_pdb $fout_options/config_t1.pdb
cp $path_pdb $fout_options/config_t2.pdb
cp $path_pdb $fout_options/config_t3.pdb
cp $path_pdb $fout_options/config_t4.pdb
python3 volgrids smiffer prot $fout_options/config_t0.pdb --config $tmp_config_options_t0
python3 volgrids smiffer prot $fout_options/config_t1.pdb --config $tmp_config_options_t1
python3 volgrids smiffer prot $fout_options/config_t2.pdb --config $tmp_config_options_t2
python3 volgrids smiffer prot $fout_options/config_t3.pdb --config $tmp_config_options_t3
python3 volgrids smiffer prot $fout_options/config_t4.pdb --config $tmp_config_options_t4

##### bonus: pocket sphere run
python3 volgrids smiffer prot $fpdb/1iqj.pdb -o $fout_pocketsphere -s 4.682 21.475 7.161 14.675 --config $tmp_config_options_t3
cp $fpdb_orig/1iqj.pdb $fout_pocketsphere/


rm -f $tmp_config_options_t0 $tmp_config_options_t1 $tmp_config_options_t2 $tmp_config_options_t3 $tmp_config_options_t4
