# Volumetric Potentials Calculator
This is the alpha version of the SMIFFER package. New features are first tested here before being implemented in the main package.

<!-- ----------------------------------------------------------------------- -->
## Setup
### Requirements
- Numpy
- Scipy
- MDAnalysis
- h5py

### Installing with conda
```
conda env create -f environment.yml
conda activate volgrids
```


<!-- ----------------------------------------------------------------------- -->
## Usage
### SMIFFER calculator
Run `python3 smiffer.py [mode] [path_input] [options...]` and provide the parameters of the calculation via arguments:
  - replace `[mode]` with `prot` or `rna` according to the structure of interest.
  - replace `[path_input]` with the path to the structure file (e.g. PDB). Mandatory positional argument.
  - Optionally, replace `[options...]` with any combination of the following:
    - `-o [path_out]` where `[path_out]` is the folder where the output SMIFs should be stored. if not provided, the parent folder of the input file will be used.
    - `-t [path_traj]`  where `[path_traj]` is the path to a trajectory file (e.g. XTC) supported by MDAnalysis. This activates "traj" mode, where SMIFs are calculated for all the frames of the trajectory and saved in a CMAP-series file.",
    - `-a [path_apbs]` where `[path_apbs]` is the path to the output of APBS. An *OpenDX* file is expected. This grid will be interpolated into the shape of the other grids.
    - `-rxyz [r] [x] [y] [z]` where `[r]`, `[x]`, `[y]` and `[z]` are the float values for the radius and X,Y,Z coordinates of a sphere in space, respectively. This activates "pocket sphere" mode, where the SMIFs will only be calculated inside the sphere provided.

### Volgrid Tools

- TODO


<!-- ----------------------------------------------------------------------- -->
## Commands examples
### APBS
```
pdb2pqr --ff=AMBER data/pdb/1iqj.pdb data/pqr/1iqj.pqr --apbs-input data/apbs/1iqj.in
apbs data/apbs/1iqj.in
```

### SMIFFER
- Perform a simple pocket sphere mode (`-rxyz`) in a protein system (`prot`).
```
python3 smiffer.py prot data/pdb/1iqj.pdb -o data/smifs -rxyz 14.675 4.682 21.475 7.161
```

- Perform a whole mode in an RNA system (`rna`) considering APBS data (`-a`).
```
python3 smiffer.py rna data/pdb/5bjo.pdb -a data/apbs/5bjo.pqr.dx -o data/smifs
```

- Perform a trajectory mode (`-t`) in an RNA system (`rna`). Note that this is only implemented for whole mode at the moment.
```
python3 smiffer.py rna data/tests/04-traj/7vki.pdb -t data/tests/04-traj/7vki.xtc -o data/tests/04-traj
```

### Volgrid Tools
- TODO


<!-- ----------------------------------------------------------------------- -->
## Benchmark
A benchmark of 10 protein-ligand and 10 rna-ligand complexes is provided at [this location](https://drive.google.com/file/d/1o1jR4RhXlIL0Jg3m0twrpbiTV7eIGZ38/view?usp=sharing), in the form of PDB and PQR input files.


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
Follow this instructions to visualize the atomic and SMIF trajectories simultaneously in Chimera. ChimeraX is recommended.
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
1) Open the PDB and load the atom trajectory into it (in ChimeraX, simply drag the files into the window).
2) Open the CMAP file in a similar way.
3) (Optional) Open the `smooth_md.py` script.
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


<!-- ----------------------------------------------------------------------- -->
## TODO
* fix concept behind hbonds that shouldn't be allowed to rotate
  * protein_backbone
  * TRP PRO HIS ARG ASN GLN
  * U C A G
* check kernel bug when running tiny systems
* making the format for the molecule params
* check what happens if performing "fix-cmap" operation when cmap input and output are the same file
* implement the fixing operation directy on "packing", to ensure that packed frames have the same resolution (add flag to override this behavior)
* implement: raise an error if a format file is opened with the wrong function
