# Volumetric Grids (VolGrids)
VolGrids is a framework for volumetric calculations, with emphasis in biological molecular systems. Some tools are also provided:
  - **APBS** via `volgrids apbs`. Requires installing [APBS](#installation-ubuntu).
  - **SMIF Calculator** via `volgrids smiffer`
  - **Volumetric Energy INSpector (VEINS)** via `volgrids veins`. WORK IN PROGRESS
  - **Volgrid Tools** via `volgrids vgtools`.

You can read more in their respective sections.

## QuickStart
```bash
pip install volgrids
volgrids --help
```

### Without installing the package
```bash
git clone https://github.com/DiegoBarMor/volgrids
cd volgrids
pip install -r requirements.txt
python3 volgrids --help
```


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- SETUP -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Setup
## Requirements
- Grid operations are based on **NumPy** arrays.
- Molecular structures and trajectories data are parsed by **MDAnalysis**.
- CMAP files are parsed by **h5py**.


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
### APBS (optional)
Follow the instructions from [here](#installation-ubuntu).


<!-- ----------------------------------------------------------------------- -->
## Usage
### Running the CLI utilities (without VolGrids installed)
You can use the tools provided by VolGrids without installing it, by calling any of the scripts in the root directory of this repository (it doesn't have to be the current directory, you can call them from anywhere). Leave `[options...]` empty to read more about the available options.

```bash
python3 volgrids apbs    [options...]
python3 volgrids smiffer [options...]
python3 volgrids veins   [options...]
python3 volgrids vgtools [options...]
```

<!-- ----------------------------------------------------------------------- -->
### Running the CLI utilities (VolGrids installed)
- VolGrids can be installed as a package via pip:
```bash
pip install volgrids
```

- (Or alternatively from its repository):
```bash
# your current directory must be the root directory of volgrids repository
pip install .
rm -rf build volgrids.egg-info # optional cleanup
```

- Then, it can be run from anywhere via:
```bash
volgrids apbs    [options...]
volgrids smiffer [options...]
volgrids veins   [options...]
volgrids vgtools [options...]
```


<!-- ----------------------------------------------------------------------- -->
### Running the tests
Follow the instructions at the [test data repo](https://github.com/DiegoBarMor/volgrids-testdata).


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- ------------------------------- SMIFFER ------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Statistical Molecular Interaction Fields (SMIF) Calculator
This is an implementation of the [Statistical Molecular Interaction Fields (SMIF)](https://www.biorxiv.org/content/10.1101/2025.04.16.649117v1) method.

## Usage
Run `python3 volgrids smiffer [mode] [path_structure] [options...]` and provide the parameters of the calculation via arguments:
  - replace `[mode]` with `prot`, `rna` or `ligand` according to the structure of interest.
  - replace `[path_structure]` with the path to the structure file (e.g. PDB). Mandatory positional argument.
  - Optionally, replace `[options...]` with any combination of the following:
    - `-o [folder_out]` where `[folder_out]` is the folder where the output SMIFs should be stored. if not provided, the parent folder of the input file will be used.
    - `-t [path_traj]`  where `[path_traj]` is the path to a trajectory file (e.g. XTC) supported by MDAnalysis. This activates "traj" mode, where SMIFs are calculated for all the frames of the trajectory and saved in a CMAP-series file.
    - `-a (path_apbs)` where `(path_apbs)` is the path to a cached output of APBS. Prevents the automatic calculation of APBS performed by volgrids' smiffer.
    - `-s [x] [y] [z] [r]` where `[x]`, `[y]`, `[z]` and `[r]` are the float values for the X,Y,Z coordinates and the radius of a sphere in space, respectively. This activates "pocket sphere" mode, where the SMIFs will only be calculated inside the sphere provided.
    - `-b [path_table]` where `[path_table]` is the path to a *.chem* table file to use for ligand mode, or to override the default macromolecules' tables. This flag is mandatory for "ligand" mode.
    - `-c [path_config]` where `[path_config]` is the path to a configuration file with global settings, to override the default settings (e.g. `config_volgrids.ini`).


<!-- ----------------------------------------------------------------------- -->
## Commands examples
- Calculate SMIFs for a protein system (`prot`) considering only the space inside a pocket sphere (`-s`).
```bash
python3 volgrids smiffer prot testdata/smiffer/pdb-nosolv/1iqj.pdb -s 4.682 21.475 7.161 14.675
```

- Calculate SMIFs for a whole RNA system (`rna`) considering APBS data (`-a`).
```bash
python3 volgrids smiffer rna testdata/smiffer/pdb-nosolv/5bjo.pdb -a testdata/smiffer/apbs/5bjo.pdb.mrc
```

- Calculate SMIFs for an RNA system (`rna`) along a trajectory (`-t`). Note that for "pocket sphere" mode, the same coordinates/radius are used for the whole trajectory.
```bash
python3 volgrids smiffer rna testdata/smiffer/traj/7vki.pdb -t testdata/smiffer/traj/7vki.xtc
```


<!-- ----------------------------------------------------------------------- -->
## Visualization
### Color standard
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
## APBS reference
Sample commands to obtain the electrostatic grids from [pdb2pqr](https://pdb2pqr.readthedocs.io/en/latest/) and [APBS](https://apbs.readthedocs.io/en/latest/)

### Installation (Ubuntu)
```bash
pip install pdb2pqr
# sudo apt install pdb2pqr # alternative
sudo apt-get install apbs
```

### Commands examples
- Running APBS with Volgrids (recommended).
```bash
python3 volgrids apbs testdata/smiffer/pdb-nosolv/1iqj.pdb --mrc --verbose
```

- Alternative (calling directly `apbs.sh`).
```bash
volgrids/utils/apbs.sh apbs testdata/smiffer/pdb-nosolv/1iqj.pdb --mrc --verbose
```

- Running APBS without Volgrids.
```bash
pdb2pqr --ff=AMBER testdata/smiffer/pdb-nosolv/1iqj.pdb testdata/smiffer/pqr/1iqj.pqr --apbs-input testdata/smiffer/1iqj.in
apbs testdata/smiffer/1iqj.in
```



<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- VEINS -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Volumetric Energy INSpector (VEINS)
This tool allows to visualize interaction energies in space by portraying them as a volumetric grid. Apart from the usual structure/trajectory files (PDB, XTC...), a CSV with energy values and the indices of the atoms/residues involved must be given. Interactions between 2, 3 and 4 particles are supported and represented accordingly

## Usage
Run `python3 volgrids veins [mode] [path_structure] [path_csv] [options...]` and provide the parameters of the calculation via arguments:
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
Run `python3 volgrids vgtools [mode] [options...]` and provide the parameters of the calculation via arguments.
  - Replace `[mode]` with one of the following available modes:
    - `convert`: Convert grid files between formats.
    - `pack`: Pack multiple grid files into a single CMAP series-file.
    - `unpack`: Unpack a CMAP series-file into multiple grid files.
    - `average`: Average all grids in a CMAP series-file into a single grid.
    - `fix_cmap`: Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.
    - `compare`: Compare two grid files by printing the number of differing points and their accumulated difference.
  - `[options...]` will depend on the mode, check the respective help string for more information (run `python3 volgrids vgtools [mode] -h`).



<!-- ----------------------------------------------------------------------- -->
