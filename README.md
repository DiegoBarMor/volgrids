# Volumetric Grids
This is a framework for volumetric calculations, with emphasis in biological molecular systems. Two tools are also provided: **SMIF Calculator** (`smiffer.py`) and **Volgrid Tools** (`vgtools.py`). You can read more in their respective sections.


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- -------------------------------- SETUP -------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Setup
## Requirements
This framework only has 2 Python dependencies:
- **MDAnalysis** for parsing structure and trajectories data. Installing this should also install **NumPy**, also needed by volgrids.
- **h5py** for parsing CMAP files.


<!-- ----------------------------------------------------------------------- -->
## Option 1: Setting up a Conda environment
### Automatic
```
conda env create -f environment/conda.yml
conda activate volgrids
```

### Manual
```
conda create --name volgrids -y
conda activate volgrids
conda install python -y
conda install -c conda-forge mdanalysis -y
```


<!-- ----------------------------------------------------------------------- -->
## Option 2: Simple setup with PIP
### Automatic
```
pip install -r environment/requirements.txt
```

### Manual
```
pip install mdanalysis h5py
```


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- ------------------------------- SMIFFER ------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# SMIF Calculator
A custom implementation of the [Statistical Molecular Interaction Fields (SMIF)](https://www.biorxiv.org/content/10.1101/2025.04.16.649117v1) method is provided as `smiffer.py`.

## Usage
Run `python3 smiffer.py [mode] [path_structure] [options...]` and provide the parameters of the calculation via arguments:
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
pdb2pqr --ff=AMBER testdata/_input/pdb/1iqj.pdb testdata/_input/pqr/1iqj.pqr --apbs-input testdata/_input/1iqj.in
apbs testdata/_input/1iqj.in
```

- Calculate SMIFs for a protein system (`prot`) considering only the space inside a pocket sphere (`-rxyz`).
```
python3 smiffer.py prot testdata/_input/pdb/1iqj.pdb -rxyz 14.675 4.682 21.475 7.161
```

- Calculate SMIFs for a whole RNA system (`rna`) considering APBS data (`-a`).
```
python3 smiffer.py rna testdata/_input/pdb/5bjo.pdb -a testdata/_input/apbs/5bjo.pqr.dx
```

- Calculate SMIFs for an RNA system (`rna`) along a trajectory (`-t`). Note that this is only implemented for "whole" mode at the moment.
```
python3 smiffer.py rna testdata/smiffer/traj/7vki.pdb -t testdata/smiffer/traj/7vki.xtc
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
| HB Acceptors    | Orange     | 1,0.5,0    | 255,128,0  | FF8000 |
| HB Donors       | Violet     | 0.7,0,1    | 179,0,255  | B300FF |

### Visualizing CMAP trajectories in Chimera
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


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- ---------------------------- VOLGRID TOOLS ---------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# Volgrid Tools
Collection of utilities for manipulating DX, MRC and CMAP grids.

## Usage
Run `python3 vgtools.py [mode] [options...]` and provide the parameters of the calculation via arguments.
  - Replace `[mode]` with `convert`, `pack`, `unpack` or `fix-cmap`. Available modes:
    - `convert`: Convert grid files between formats.
    - `pack`: Pack multiple grid files into a single CMAP series-file.
    - `unpack`: Unpack a CMAP series-file into multiple grid files.
    - `fix-cmap`: Ensure that all grids in a CMAP series-file have the same resolution, interpolating them if necessary.
  - `[options...]` will depend on the mode, check the respective help string for more information (run `python3 vgtools.py [mode]` with no more arguments).


<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- ----------------------------------------------------------------------- -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
# TODO
* fix concept behind hbonds that shouldn't be allowed to rotate
  * protein_backbone
  * TRP PRO HIS ARG ASN GLN
  * U C A G
* check what happens if performing "fix-cmap" operation when cmap input and output are the same file
* implement the fixing operation directy on "packing", to ensure that packed frames have the same resolution (add flag to override this behavior)
* implement: raise an error if a format file is opened with the wrong function
* improve IniParser.parse_str
* provide download links for the testdata folder
* add annotations, docstrings and overall cleaning


<!-- ----------------------------------------------------------------------- -->
