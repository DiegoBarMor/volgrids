# Volumetric Grids (VolGrids)
VolGrids is a framework for volumetric calculations, with emphasis in biological molecular systems. The following applications are provided:
  - [**SMIF Calculator**](#statistical-molecular-interaction-fields-smif-calculator) via `volgrids smiffer`. This is an implementation of the [Statistical Molecular Interaction Fields (SMIF)](https://www.biorxiv.org/content/10.1101/2025.04.16.649117v1) method.
  - [**Smiffer Utilities**](#smiffer-utilities)  via `volgrids smutils`. Utilities related to more advanced SMIF usage.
  - **APBS** via `volgrids apbs`. Requires installing [APBS](#installation-ubuntu).
  - [**Volgrid Tools**](#volgrid-tools) via `volgrids vgtools`. Collection of utilities for manipulating DX, MRC, CCP4 and CMAP grids.

You can read more in the [subcommands summary](#summary-of-subcommands).

## QuickStart
```bash
pip install volgrids
volgrids --help
```

### Without installing the package
```bash
git clone --depth 1 https://github.com/DiegoBarMor/volgrids
cd volgrids
pip install -r requirements.txt
python3 volgrids --help # required vendors will automatically be downloaded
```


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- SETUP -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Setup
## Requirements
### Hard requirements
- Grid operations are based on **NumPy** arrays.
- Molecular structures and trajectories data are parsed by [**MDAnalysis**](https://github.com/MDAnalysis/mdanalysis).
- CMAP files are parsed by **h5py**.

### Optional requirements
- [**APBS**]: Follow the instructions from [here](#installation-ubuntu).
- [**rnapolis**](https://github.com/tzok/rnapolis-py) (`pip install rnapolis`) for running residue-selection utilities from `smutils`.

### Vendors
The following are dependencies that are included as vendors in `volgrid`'s pip distribution. If you don't want to use pip, it will be fetched automatically the first time you run `volgrids`. Alternatively, run `bash scripts/_prepare.sh` in your local copy of the `volgrids` repo to fetch the vendor packages manually.
- [**freyacli**](https://github.com/DiegoBarMor/freyacli) for CLI management
- [**molutils**](https://github.com/DiegoBarMor/molutils) for utilities dealing with PDB files.


<!-- ----------------------------------------------------------------------- -->
### Option 1: Setting up a Conda environment
#### Automatic
```bash
conda env create -f environment.yml
conda activate volgrids
```

#### Manual
```bash
conda create --name volgrids -y
conda activate volgrids
conda install python -y
conda install -c conda-forge mdanalysis h5py -y
```


<!-- ----------------------------------------------------------------------- -->
### Option 2: Simple setup with PIP
#### Automatic
```bash
pip install -r requirements.txt
```

#### Manual
```bash
pip install mdanalysis h5py
```


<!-- ----------------------------------------------------------------------- -->
## Running the CLI utilities
### Without VolGrids installed
You can use the tools provided by VolGrids without installing it, by calling any of the scripts in the root directory of this repository (it doesn't have to be the current directory, you can call them from anywhere). Leave `[options...]` empty to read more about the available options.

```bash
python3 volgrids apbs    [options...]
python3 volgrids smiffer [options...]
python3 volgrids smutils [options...]
python3 volgrids vgtools [options...]
```

<!-- ----------------------------------------------------------------------- -->
### VolGrids installed
- VolGrids can be installed as a package via pip:
```bash
pip install volgrids
```

- (Or alternatively from its repository):
```bash
# your current directory must be the root directory of volgrids repository
git clone --depth 1 https://github.com/DiegoBarMor/volgrids
cd volgrids
pip install .
rm -rf build volgrids.egg-info # optional cleanup
```

- Then, it can be run from anywhere via:
```bash
volgrids apbs    [options...]
volgrids smiffer [options...]
volgrids smutils [options...]
volgrids vgtools [options...]
```


<!-- ----------------------------------------------------------------------- -->
## Running the tests
Follow the instructions at the [test data repo](https://github.com/DiegoBarMor/volgrids-testdata).


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- USAGE -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Usage
## Summary of subcommands
Run the subcommands without any further arguments to read more about their specific usage.
- `volgrids apbs`: Generate raw APBS potential grids for biomolecular structures.
- `volgrids smiffer prot`: Calculate SMIFs for protein structures.
- `volgrids smiffer rna`: Calculate SMIFs for RNA structures.
- `volgrids smiffer ligand`: Calculate SMIFs for ligand structures. A .chem table must be provided.
- `volgrids smutils resids_nobp`: Print the set of non-base-paired residue indices in a given RNA structure. A residue is considered non-base-paired if it does not form a canonical base pair (UA, CG) with any other residue. Requires rnapolis.
- `volgrids smutils resids_nostk`: Print the set of residue indices in a given RNA structure for residues that aren't participating in 2 stacking interactions. Requires rnapolis.
- `volgrids smutils histogram`: Plot the distribution histogram of non-zero voxel values in a grid file. For multi-frame CMAP files, all frames are concatenated unless --key is specified. Percentiles (p50, p75, p90, p95, p99, p99.9) are printed to the console. Requires matplotlib.
- `volgrids smutils occupancy prot`: Calculate OGs for protein structures.
- `volgrids smutils occupancy rna`: Calculate OGs for RNA structures.
- `volgrids smutils occupancy ligand`: Calculate OGs for ligand structures. A .chem table must be provided.
- `volgrids vgtools convert`: Convert grid files between formats.
- `volgrids vgtools pack`: Pack multiple grid files into a single CMAP series-file.
- `volgrids vgtools unpack`: Unpack a CMAP series-file into multiple grid files.
- `volgrids vgtools fix_cmap`: Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.
- `volgrids vgtools average`: Average all grids in a CMAP series-file into a single grid.
- `volgrids vgtools summary`: Print a summary of the grid file (format, dimensions, resolution, etc.) to the console.
- `volgrids vgtools compare`: Compare two grid files by printing the number of differing points and their accumulated difference.
- `volgrids vgtools rotate`: Rotate a grid file by 3 angles, along the xy, yz and xz planes (in degrees).
- `volgrids vgtools op abs`: Element-wise absolute value of the grid's points.
- `volgrids vgtools op add`: Element-wise addition of the grids' points.
- `volgrids vgtools op sub`: Element-wise substraction of the grids' points. Order is "grid_0 - grid_1".
- `volgrids vgtools op mul`: Element-wise multiplication of the grids' points.
- `volgrids vgtools op div`: Element-wise division of the grids' points. Order is "grid_0 / grid_1". Dubious operation, as grids have usually many 0-valued points; make sure to preprocess the denominator grid appropriately.


<!-- ----------------------------------------------------------------------- -->
## Commands examples (smiffer)
- Calculate SMIFs for a protein system (`prot`) considering only the space inside a pocket sphere (`-s`).
```bash
volgrids smiffer prot testdata/smiffer/pdb_clean/1iqj.pdb -s 4.682 21.475 7.161 14.675
```

- Calculate SMIFs for a whole RNA system (`rna`) considering APBS data (`-a`).
```bash
volgrids smiffer rna testdata/smiffer/pdb_clean/5bjo.pdb -a testdata/smiffer/apbs/5bjo.pdb.mrc
```

- Calculate SMIFs for an RNA system (`rna`) along a trajectory (`-t`). Note that for "pocket sphere" mode, the same coordinates/radius are used for the whole trajectory.
```bash
volgrids smiffer rna testdata/smiffer/traj/7vki.pdb -t testdata/smiffer/traj/7vki.xtc
```

- Combine `resids_nonpb` with smiffer's `--resids` (`-i`) flag to have more polished hbond results, e.g.:
```bash
python3 volgrids smiffer rna 1akx.pdb \
    -i "$(python3 volgrids smutils resids_nobp 1akx.pdb)" \
    -c DO_SMIF_STACKING=false DO_SMIF_HYDROPHOBIC=false DO_SMIF_HYDROPHILIC=false \
        DO_SMIF_HBA=true DO_SMIF_HBD=true HBONDS_ONLY_NUCLEOBASE=true
```


<!-- ----------------------------------------------------------------------- -->
## Visualization
### Color standard (smiffer)
| Potential       | Color      | RGB 0-1    | RGB 0-255  | HEX    |
|-----------------|------------|------------|------------|--------|
| APBS -          | Red        | 1,0,0      | 255,0,0    | FF0000 |
| APBS +          | Blue       | 0,0,1      | 0,0,255    | 0000FF |
| HB Acceptors    | Violet     | 0.7,0,1    | 179,0,255  | B300FF |
| HB Donors       | Orange     | 1,0.5,0    | 255,128,0  | FF8000 |
| Hydrophilic (-) | Light Blue | 0.3,0.85,1 | 77,217,255 | 4DD9FF |
| Hydrophobic (+) | Yellow     | 1,1,0      | 255,255,0  | FFFF00 |
| Stacking        | Green      | 0,1,0      | 0,255,0    | 00FF00 |


### MRC/CCP4 data in ChimeraX
Use this command when visualizing MRC/CCP4 data with negative values in ChimeraX (replace `1` with the actual number of the model).
```
volume #1 capFaces false
```


### CMAP trajectories in ChimeraX
Follow these instructions to visualize the atomic and SMIF trajectories simultaneously in ChimeraX.
1) Open the PDB and load the atom trajectory into it (in ChimeraX, simply drag the files into the window).
2) Open the CMAP file in a similar way.
3) Start the playback by using this ChimeraX command. The numbers specified would change if dealing with multiple structures/cmaps. Examples:
```
coordset #1; vseries play #2
coordset #1 pauseFrames 5; vseries play #2 pauseFrames 5
coordset #1 pauseFrames 5; vseries play #2 pauseFrames 5; vseries play #3 pauseFrames 5
```
4) Use this ChimeraX command to stop the playback. The ids used must match the previous command.
```
coordset stop #1; vseries stop #2
```


#### Smooth Trajectories
1) Load the PDB and the trajectory files into it ChimeraX (e.g. drag the files into the window).
2) Load the CMAP file in a similar way.
3) (Optional) Load the `smooth_md.py` script (again, can be done by dragging it into ChimeraX).
4) Start the playback by using this ChimeraX command. The numbers specified would change if dealing with multiple structures/cmaps. Examples:
```
coordset #1 pauseFrames 10; vop morph #2 playStep 0.0005 frames 2000 modelId 3
coordset #1 pauseFrames 20; vop morph #2 playStep 0.00025 frames 4000 modelId 3
coordset #1 pauseFrames 20; vop morph #2 playStep 0.00025 frames 4000 modelId 4; vop morph #3 playStep 0.00025 frames 4000 modelId 5
```
4) Use this ChimeraX command to stop the playback. The ids used must match the previous command.
```
coordset stop #1; vseries stop #2
```
Note that this time, the morph can be paused manually with the slider button (is there a command equivalent?)

#### Other useful ChimeraX commands
```
volume level 0.5
volume transparency 0.5
volume showOutlineBox true

save trajectory.pdb models #1 allCoordsets true
cartoon style width 0.2 thickness 0.1
size stickradius 0.1
```


<!-- ----------------------------------------------------------------------- -->
# APBS reference
Sample commands to obtain the electrostatic grids from [pdb2pqr](https://pdb2pqr.readthedocs.io/en/latest/) and [APBS](https://apbs.readthedocs.io/en/latest/)

## Installation (Ubuntu)
```bash
pip install pdb2pqr
# sudo apt install pdb2pqr # alternative
sudo apt-get install apbs
```

## Commands examples
- Running APBS with Volgrids (recommended).
```bash
volgrids apbs testdata/smiffer/pdb_clean/1iqj.pdb --mrc --verbose
```

- Alternative (calling directly `apbs.sh`).
```bash
bash volgrids/apbs/apbs.sh testdata/smiffer/pdb_clean/1iqj.pdb --mrc --verbose
```

- Running APBS without Volgrids.
```bash
pdb2pqr --ff=AMBER testdata/smiffer/pdb_clean/1iqj.pdb testdata/smiffer/pqr/1iqj.pqr --apbs-input testdata/smiffer/1iqj.in
apbs testdata/smiffer/1iqj.in
```

- Reference for the `.in` file: https://apbs.readthedocs.io/en/latest/using/input/old/elec/mg-auto.html


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
