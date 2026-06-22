#!/bin/bash
set -eu

echo
echo ">>> TEST SMIFFER 6: CavityFinder"

fpdb="testdata/smiffer/pdb_clean"
fpdb_orig="testdata/smiffer/pdb_orig"
fout="testdata/smiffer/cavities"
fout_benchmark=$fout/benchmark
fout_pocketsphere=$fout/pocket_sphere
fout_options=$fout/options
rm -rf $fout
mkdir -p $fout_benchmark $fout_pocketsphere $fout_options

conf_cavities() { # trim_occ, trim_cav, save_trim, save_cav, threshold
    echo "GRID_FORMAT_OUTPUT=CMAP DO_TRIMMING_OCCUPANCY=$1 DO_TRIMMING_CAVITIES=$2 SAVE_TRIMMING_MASK=$3 SAVE_CAVITIES=$4 TRIMMING_CAVITIES_THRESHOLD=$5"
}
conf_no_smifs="GRID_FORMAT_OUTPUT=CMAP DO_SMIF_STACKING=False DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"
conf_just_stak="GRID_FORMAT_OUTPUT=CMAP DO_SMIF_STACKING=True DO_SMIF_HBA=False DO_SMIF_HBD=False DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"
conf_simple_smifs="GRID_FORMAT_OUTPUT=CMAP DO_SMIF_STACKING=True DO_SMIF_HBA=True DO_SMIF_HBD=True DO_SMIF_HYDROPHOBIC=False DO_SMIF_HYDROPHILIC=False DO_SMIF_APBS=False"


############################# BENCHMARK SYSTEMS
conf_benchmark="$conf_no_smifs $(conf_cavities True False False True 3)"

run_benchmarks() {
    local names=$1
    for name in $names; do
        for i in 1 2 3; do
            python3 volgrids smiffer "$fpdb/$name.pdb" -o $fout_benchmark --config "$conf_benchmark" CAVITIES_NPASSES=$i
            mv "$fout_benchmark/$name.cav.cmap" $fout_benchmark/npasses_$i.cmap
        done

        python3 volgrids smiffer "$fpdb/$name.pdb" -o $fout_benchmark --config "$conf_simple_smifs" \
            CAVITIES_WEIGHT=1.0 GRID_FORMAT_OUTPUT=MRC CAVITIES_NPASSES=2

        for smif in stk hba hbd; do
            mv "$fout_benchmark/$name.$smif.mrc" "$fout_benchmark/$name.$smif.weighted.mrc"
        done

        python3 volgrids smiffer "$fpdb/$name.pdb" -o $fout_benchmark --config "$conf_simple_smifs" SAVE_TRIMMING_MASK=true
        (
            shopt -s nullglob
            cmaps=( "$fout_benchmark"/npasses*.cmap "$fout_benchmark/$name".*.cmap )
            ### no need to include the last CMAP created, as "pack" will append the other files to this one
            python3 volgrids vgtools pack "${cmaps[@]}" -o "$fout_benchmark/$name.cmap"
        )
        rm -f "$fout_benchmark/$name".*.weighted.mrc

        cp "$fpdb_orig/$name.pdb" $fout_benchmark/

    done
    rm -f $fout_benchmark/npasses_*
}

run_benchmarks "1bg0 1eby" # 1ehe 1h7l 1iqj 1ofz 3dd0 3ee4 5m9w 6e9a"   # prot
run_benchmarks "1akx 1i9v" # 2esj 4f8u 5bjo 5kx9 6tf3 7oax0 7oax1 8eyv" # rna


############################# OPTIONS
path_pdb="$fpdb/1iqj.pdb"

######################### booleans
conf_options_00="$conf_just_stak $(conf_cavities False True  True True  3)" # DO_TRIMMING_OCCUPANCY=false, so cavities should be ignored
conf_options_01="$conf_just_stak $(conf_cavities True  False True False 3)" # cavities shouldn't be saved and shouldn't affect the trimming
conf_options_02="$conf_just_stak $(conf_cavities True  True  True False 3)" # cavities shouldn't be saved but should affect the trimming
conf_options_03="$conf_just_stak $(conf_cavities True  False True True  3)" # cavities should be saved but shouldn't affect the trimming
conf_options_04="$conf_just_stak $(conf_cavities True  True  True True  3)" # cavities should be saved and should affect the trimming

cp $path_pdb $fout_options/config_00.pdb
cp $path_pdb $fout_options/config_01.pdb
cp $path_pdb $fout_options/config_02.pdb
cp $path_pdb $fout_options/config_03.pdb
cp $path_pdb $fout_options/config_04.pdb
python3 volgrids smiffer $fout_options/config_00.pdb --config "$conf_options_00"
python3 volgrids smiffer $fout_options/config_01.pdb --config "$conf_options_01"
python3 volgrids smiffer $fout_options/config_02.pdb --config "$conf_options_02"
python3 volgrids smiffer $fout_options/config_03.pdb --config "$conf_options_03"
python3 volgrids smiffer $fout_options/config_04.pdb --config "$conf_options_04"


######################### threshold
conf_options_t0="$conf_just_stak $(conf_cavities True True True True 0)"
conf_options_t1="$conf_just_stak $(conf_cavities True True True True 1)"
conf_options_t2="$conf_just_stak $(conf_cavities True True True True 2)"
conf_options_t3="$conf_just_stak $(conf_cavities True True True True 3)"
conf_options_t4="$conf_just_stak $(conf_cavities True True True True 4)"

cp $path_pdb $fout_options/config_t0.pdb
cp $path_pdb $fout_options/config_t1.pdb
cp $path_pdb $fout_options/config_t2.pdb
cp $path_pdb $fout_options/config_t3.pdb
cp $path_pdb $fout_options/config_t4.pdb
python3 volgrids smiffer $fout_options/config_t0.pdb --config "$conf_options_t0"
python3 volgrids smiffer $fout_options/config_t1.pdb --config "$conf_options_t1"
python3 volgrids smiffer $fout_options/config_t2.pdb --config "$conf_options_t2"
python3 volgrids smiffer $fout_options/config_t3.pdb --config "$conf_options_t3"
python3 volgrids smiffer $fout_options/config_t4.pdb --config "$conf_options_t4"

##### bonus: pocket sphere run
python3 volgrids smiffer $fpdb/1iqj.pdb -o $fout_pocketsphere -s 4.682 21.475 7.161 14.675 --config "$conf_options_t3"
cp $fpdb_orig/1iqj.pdb $fout_pocketsphere/
