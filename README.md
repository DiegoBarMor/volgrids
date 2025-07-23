# Volumetric Grids (VolGrids)
This is a framework for volumetric calculations, with emphasis in biological molecular systems. Three tools are also provided: **SMIF Calculator** (`run/smiffer.py`), **Volumetric Energy INSpector (VEINS)** (`run/veins.py`) and **Volgrid Tools** (`run/vgtools.py`). You can read more in their respective sections.


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- SETUP -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Setup
## Requirements
This framework only has 2 Python dependencies:
- **MDAnalysis** for parsing structure and trajectories data. Installing this should also install **NumPy**, also needed by volgrids.
- **h5py** for parsing CMAP files.


<!-- ----------------------------------------------------------------------- -->
### Option 1: Setting up a Conda environment
#### Automatic
```
conda env create -f environment/conda.yml
conda activate volgrids
```

#### Manual
```
conda create --name volgrids -y
conda activate volgrids
conda install python -y
conda install -c conda-forge mdanalysis -y
```


<!-- ----------------------------------------------------------------------- -->
### Option 2: Simple setup with PIP
#### Automatic
```
pip install -r environment/requirements.txt
```

#### Manual
```
pip install mdanalysis h5py
```


<!-- ----------------------------------------------------------------------- -->
## Usage
### Without installing VolGrids
You can use the tools provided by VolGrids without installing it, by calling any of the scripts in the root directory of this repository (it doesn't have to be the current directory, you can call them from anywhere). Leave `[options...]` empty to read more about the available options.

- **SMIF Calculator:**
```
python3 run/smiffer.py [options...]
```

- **Volumetric Energy INSpector (VEINS):**
```
run/python3 veins.py [options...]
```

- **Volgrid Tools:**
```
python3 run/vgtools.py [options...]
```


<!-- ----------------------------------------------------------------------- -->
### Using VolGrids as a package
You can install VolGrids as a package and import it from your own scripts. For installing with pip:
```
# your current directory should be the root folder of this repository
pip install .
rm -rf build volgrids.egg-info # optional cleanup
```


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- ------------------------------- SMIFFER ------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Statistical Molecular Interaction Fields (SMIF) Calculator
This is a custom implementation of the [Statistical Molecular Interaction Fields (SMIF)](https://www.biorxiv.org/content/10.1101/2025.04.16.649117v1) method.

## Usage
Run `python3 run/smiffer.py [mode] [path_structure] [options...]` and provide the parameters of the calculation via arguments:
  - replace `[mode]` with `prot`, `rna` or `ligand` according to the structure of interest.
  - replace `[path_structure]` with the path to the structure file (e.g. PDB). Mandatory positional argument.
  - Optionally, replace `[options...]` with any combination of the following:
    - `-o [folder_out]` where `[folder_out]` is the folder where the output SMIFs should be stored. if not provided, the parent folder of the input file will be used.
    - `-t [path_traj]`  where `[path_traj]` is the path to a trajectory file (e.g. XTC) supported by MDAnalysis. This activates "traj" mode, where SMIFs are calculated for all the frames of the trajectory and saved in a CMAP-series file.
    - `-a [path_apbs]` where `[path_apbs]` is the path to the output of APBS. An *OpenDX* file is expected. This grid will be interpolated into the shape of the other grids.
    - `-rxyz [r] [x] [y] [z]` where `[r]`, `[x]`, `[y]` and `[z]` are the float values for the radius and X,Y,Z coordinates of a sphere in space, respectively. This activates "pocket sphere" mode, where the SMIFs will only be calculated inside the sphere provided.
    - `-b [path_table]` where `[path_table]` is the path to a *.chem* table file to use for ligand mode, or to override the default macromolecules' tables. This flag is mandatory for "ligand" mode.
    - `-c [path_config]` where `[path_config]` is the path to a configuration file with global settings, to override the default settings from `config.ini`.


<!-- ----------------------------------------------------------------------- -->
## Commands examples
- Sample commands to obtain the electrostatic grids from [pdb2pqr](https://pdb2pqr.readthedocs.io/en/latest/) and [APBS](https://apbs.readthedocs.io/en/latest/)
```
pdb2pqr --ff=AMBER testdata/_input/pdb-nosolv/1iqj.pdb testdata/_input/pqr/1iqj.pqr --apbs-input testdata/_input/1iqj.in
apbs testdata/_input/1iqj.in
```

- Calculate SMIFs for a protein system (`prot`) considering only the space inside a pocket sphere (`-rxyz`).
```
python3 run/smiffer.py prot testdata/_input/pdb-nosolv/1iqj.pdb -rxyz 14.675 4.682 21.475 7.161
```

- Calculate SMIFs for a whole RNA system (`rna`) considering APBS data (`-a`).
```
python3 run/smiffer.py rna testdata/_input/pdb-nosolv/5bjo.pdb -a testdata/_input/apbs/5bjo.pqr.dx
```

- Calculate SMIFs for an RNA system (`rna`) along a trajectory (`-t`). Note that this is only implemented for "whole" mode at the moment.
```
python3 run/smiffer.py rna testdata/smiffer/traj/7vki.pdb -t testdata/smiffer/traj/7vki.xtc
```


<!-- ----------------------------------------------------------------------- -->
## Benchmark
A benchmark of 10 protein-ligand and 10 rna-ligand complexes is provided at [this location](https://drive.google.com/file/d/1o1jR4RhXlIL0Jg3m0twrpbiTV7eIGZ38/view?usp=sharing), in the form of PDB and PQR input files.
<!-- TODO: update this with a new link to the testdata folder -->


<!-- ----------------------------------------------------------------------- -->
## Visualization
### Color standard
| Potential       | Color      | RGB 0-1    | RGB 0-255  | HEX    |
|-----------------|------------|------------|------------|--------|
| APBS +          | Blue       | 0,0,1      | 0,0,255    | 0000FF |
| APBS -          | Red        | 1,0,0      | 255,0,0    | FF0000 |
| Hydrophobic (+) | Yellow     | 1,1,0      | 255,255,0  | FFFF00 |
| Hydrophilic (-) | Light Blue | 0.3,0.85,1 | 77,217,255 | 4DD9FF |
| HB Donors       | Orange     | 1,0.5,0    | 255,128,0  | FF8000 |
| HB Acceptors    | Violet     | 0.7,0,1    | 179,0,255  | B300FF |
| Stacking        | Green      | 0,1,0      | 0,255,0    | 00FF00 |


### MRC/CCP4 data in Chimera
Use this command when visualizing MRC/CCP4 data with negative values in Chimera (replace `1` with the actual number of the model).
```
volume #1 capFaces false
```


### CMAP trajectories in Chimera
Follow these instructions to visualize the atomic and SMIF trajectories simultaneously in Chimera. ChimeraX is recommended.
1) Open the PDB and load the atom trajectory into it (in ChimeraX, simply drag the files into the window).
2) Open the CMAP file in a similar way.
3) Start the playback by using this Chimera command. The numbers specified would change if dealing with multiple structures/cmaps. Examples:
```
coordset #1; vseries play #2
coordset #1 pauseFrames 5; vseries play #2 pauseFrames 5
coordset #1 pauseFrames 5; vseries play #2 pauseFrames 5; vseries play #3 pauseFrames 5
```
4) Use this Chimera command to stop the playback. The ids used must match the previous command.
```
coordset stop #1; vseries stop #2
```


#### Smooth Trajectories
1) Load the PDB and the trajectory files into it Chimera (in ChimeraX, simply drag the files into the window).
2) Load the CMAP file in a similar way.
3) (Optional) Load the `smooth_md.py` script (again, can be done by dragging it into ChimeraX).
4) Start the playback by using this Chimera command. The numbers specified would change if dealing with multiple structures/cmaps. Examples:
```
coordset #1 pauseFrames 10; vop morph #2 playStep 0.0005 frames 2000 modelId 3
coordset #1 pauseFrames 20; vop morph #2 playStep 0.00025 frames 4000 modelId 3
coordset #1 pauseFrames 20; vop morph #2 playStep 0.00025 frames 4000 modelId 4; vop morph #3 playStep 0.00025 frames 4000 modelId 5
```
4) Use this Chimera command to stop the playback. The ids used must match the previous command.
```
coordset stop #1; vseries stop #2
```
Note that this time, the morph can be paused manually with the slider button (is there a command equivalent?)

#### Other useful Chimera commands
```
volume level 0.5
volume transparency 0.5
```



<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- VEINS -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Volumetric Energy INSpector (VEINS)
This tool allows to visualize interaction energies in space by portraying them as a volumetric grid. Apart from the usual structure/trajectory files (PDB, XTC...), a CSV with energy values and the indices of the atoms/residues involved must be given. Interactions between 2, 3 and 4 particles are supported and represented accordingly

## Usage
Run `run/python3 veins.py [mode] [path_structure] [path_csv] [options...]` and provide the parameters of the calculation via arguments:
  - replace `[mode]` with `energies`.
  - replace `[path_structure]` with the path to the structure file (e.g. PDB). Mandatory positional argument.
  - replace `[path_csv]` with the path to the energies CSV file. Mandatory positional argument. It must contain the following rows:
    - **kind**: Name of the interaction kind. All rows with the same *kind* will be used to calculate a single grid with its name.
    - **npoints**: Number of particles involved in the interaction.
    - **idxs**: Group of 0-based indices joined by `-`. These are the indices of the particles involved in the interaction. This group must contain *npoints* indices.
    - **idxs_are_residues**: Whether the indices correspond to the molecule's residues (`true`) or atoms (`false`).
    - **energy**: Value of the interaction's energy.
  - Optionally, replace `[options...]` with any combination of the following:
    - `-o [folder_out]` where `[folder_out]` is the folder where the output SMIFs should be stored. if not provided, the parent folder of the input file will be used.
    `-c [cutoff]` where `[cutoff]` is a float number. Energies below this cutoff will be ignored. Default value: 1e-3.



<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- ---------------------------- VOLGRID TOOLS ---------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Volgrid Tools
Collection of utilities for manipulating DX, MRC, CCP4 and CMAP grids.

## Usage
Run `python3 run/vgtools.py [mode] [options...]` and provide the parameters of the calculation via arguments.
  - Replace `[mode]` with one of the following available modes:
    - `convert`: Convert grid files between formats.
    - `pack`: Pack multiple grid files into a single CMAP series-file.
    - `unpack`: Unpack a CMAP series-file into multiple grid files.
    - `fix_cmap`: Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.
    - `compare`: Compare two grid files by printing the number of differing points and their accumulated difference.
  - `[options...]` will depend on the mode, check the respective help string for more information (run `python3 run/vgtools.py [mode] -h`).



<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- --------------------------------- TODO -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# TODO
* improve help strings
* add annotations, docstrings and overall cleaning
* try out cavities-finder ideas

## VGrids
* implement: raise an error if a format file is opened with the wrong function
* make a common App class
  * add tests for parameters being directly passed to the App classes (instead of parsing the CLI arguments)
* standard ini files use ; instead of # for comments

## SMIFFER
* maybe: replace the RNDS trimming with a faster method
* provide download links for the testdata folder
* change the ligand example to one that uses both NAMES_HBACCEPTORS, NAMES_HBDONORS and NAMES_HBD_FIXED
* document the .chem tables
* check if there's a bug in the peptide bond N of the test toy system peptide_no_h

## VEINS
* implement "forces" mode
* move Grid's static fields into config.ini
* add tests

## VGTools
* check what happens if performing "fix_cmap" operation when cmap input and output are the same file
* implement the fixing operation directy on "packing", to ensure that packed frames have the same resolution (add flag to override this behavior)
* mode to describe grids
* mode to perform operations on grids: abs, sum, diff, mask...
* when editing a CMAP file (be it converting it or performing an operation on it), one should be able to specify the key of the relevant grid (instead of GridIO.read_auto arbitrarily deciding to take the first key it finds in the CMAP header)


<!-- ----------------------------------------------------------------------- -->
