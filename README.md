# Volumetric Potentials Calculator
## Required packages
- Numpy
- Scipy
- MDAnalysis

## Usage
Two modes are available for this standalone calculator: **PocketSphere (PS)** and **Whole (W)**. **PS** calculates the potentials in a spherical volume defined to be a binding pocket, while **W** calculates the potentials of all the volume surrounding the macromolecule.

### Running via "main.py"
Run `main.py` and provide the parameters of the calculation via arguments:
  - **(First positional argument)**: Path to the *PDB* structure file of interest.
  - **(Second positional argument)**: Path to the folder where the output potentials should be stored.
  - `apbs` (-a) : Path to the output of APBS. An *OpenDX* file is expected.
  - `rna` (-n) : The target macromolecule of interest is nucleic. If this flag is not provided, the target macromolecule is assumed to be proteic.
  - `whole` (-w) : Activates 'whole mode'. The potentials of all the volume surrounding the macromolecule will be calculated.
  - `default-res` (-s) : Use the default resolution values from settings.py (instead of the default delta values).
  - `radius` (-r) : Radius of the pocket sphere. Will be ignored if 'whole' mode is active.
  - `xcog` (-x) : X coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.
  - `ycog` (-y) : Y coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.
  - `zcog` (-z) : Z coordinate for the center of geometry of the pocket sphere. Will be ignored if 'whole' mode is active.

### APBS
The potential grids provided by APBS can be trimmed and interpolated into the same shape of the other grids, by using the flag `-a`. In this case, the APBS pipeline must be run first to obtain the input electrostatic potential grids. Example of using APBS:
```
pdb2pqr --ff=AMBER data/pdb/1iqj.pdb data/pqr/1iqj.pqr --apbs-input data/apbs/1iqj.in
apbs data/apbs/1iqj.in
```

### Examples
- Perform a simple pocket sphere mode (`-x`, `-y`, `-z`, `-r`) in a protein system.
```
python main.py data/pdb/1iqj.pdb data/potentials -x 4.682 -y 21.475 -z 7.161 -r 14.675
```

- Perform a whole mode (`-w`) in an RNA system (`-n`) considering APBS data (`-a`).
```
python main.py data/pdb/5bjo.pdb data/potentials-whole -w -n -a data/apbs/5bjo.pqr.dx
```

## Benchmark
A benchmark of 10 protein-ligand and 10 rna-ligand complexes is provided at [this location](https://drive.google.com/file/d/1o1jR4RhXlIL0Jg3m0twrpbiTV7eIGZ38/view?usp=sharing), in the form of PDB and PQR input files.

## Color standard
| Potential       | Color      | RGB 0-1    | RGB 0-255  | HEX    |
|-----------------|------------|------------|------------|--------|
| APBS +          | Blue       | 0,0,1      | 0,0,255    | 0000FF |
| APBS -          | Red        | 1,0,0      | 255,0,0    | FF0000 |
| Hydrophobic (+) | Yellow     | 1,1,0      | 255,255,0  | FFFF00 |
| Hydrophilic (-) | Light Blue | 0.3,0.85,1 | 77,217,255 | 4DD9FF |
| HB Acceptors    | Orange     | 1,0.5,0    | 255,128,0  | FF8000 |
| HB Donors       | Violet     | 0.7,0,1    | 179,0,255  | B300FF |
