# Changelog

## [0.14.0] - 2026-06-02
- Refactored and cleaned up `.chem` tables.
    - Added distinction between RNA and DNA residues.
    - Simplified header names.

    - Atoms of a stacking group are now to be provided separated by a `-` character.
    - For residues with multiple stacking groups, they are to be provided in the same row, separated by spaces.

    - Removed the `*` wildcard for hydrophobicity atoms.
        - Names for atoms with the same hydrophobicity value can now be more conveniently grouped using the `,` char as delimiter, e.g.: `URA: O5',P,OP1,OP2,O3'=-2.02 C1',C2',C3',C4',C5',O2',O4'=-0.13 N1,C2,N3,C4,C5,C6,O2,O4=1.13`.
    - Removed the generic hydrophobicity scale for whole residues (now it's always specificed for individual atoms of a given residue).
    - Default hydrophobicity scale for proteins now ignores the common peptide atoms (i.e. `N CA C O`).

- Moved imports of some rare/heavy dependencies to inside relevant methods.
    - This way, the heavy overhead of libraries like MDAnalysis can be avoided when performing simple operation that doesn't need them (e.g. running `volgrids -h`).


## [0.13.1] - 2026-06-01
- Added `BIN` format for grids.
    - `BIN` is a straightforward binary format that packs a grid in the following way (little-endian):
        3 unsigned ints: resolution (nx, ny, nz); 3 floats: deltas; 3 floats: origin; nx*ny*nz float32 values: flattened array in C-order
    - it's meant for easy parsing in lower level languages


## [0.13.0] - 2026-05-28
- Added new subcommand group for `smutils`: `sphere`.
    - `smutils sphere find`: Find the sphere information (x,y,z,radius) for the smallest sphere that encloses all the atoms inside a selection query.
    - `smutils sphere grid`: Create a grid with a single boolean sphere. The grid encloses completely the provided structure file.
    - These subcommands are meant to be used together (see example usage in `calc_spheres()` inside `examples/spheres/run.sh`).
    - They can also consider trajectory files (for considering a single sphere moving in time).
- Some code improvements.


## [0.12.1] - 2026-05-27
- Fixed bug where trimming was considering water and other unwanted atoms.
    - This bug was introduced in `0.12.0` when the default chem files were refactored.
- Fixed bug where HBDonor OGs were including unwanted atoms.
    - The special cases for RNAs and proteins were not being taken into account.
- Added `std_dev` subcommand to `vgtools`.
    - Similar to the `average` subcommand, it allows to obtain the standard deviation grid for multiple grids in a CMAP file.


## [0.12.0] - 2026-05-22
- CLI changes in `smiffer` (**IMPORTANT**: breaks compatibility with calls to older versions of the smiffer's command).
    - Renamed the `-b`/`--table` flag to `-e`/`--chem` to better represent its intent (i.e. use it for passing a custom `.chem` file).
    - Removed the `prot`/`rna`/`macro`/`ligand` subcommand for using smiffer.
        - Motivation: the subcommand was cumbersome and error prone.
        - Now, a single `default.chem` file contains the data for both protein and nucleic structures. By default, SMIFs will be calculated for all peptide and nucleic acids found in the provided PDB (the user must perform any relevant preprocessing if they want to focus on only one macromolecule type).
            - Ligand mode is achieved by simply using the `-e`/`--chem` flag.

- CLI changes in `vgtools` / `smutils`.
    - moved `histogram` subcommand to vgtools.
    - reorganized vgtools and smutils subcommands order.


## [0.11.2] - 2026-05-19
- Added new experimental features to `smutils`.
    - `chemgen` subcommand for automatic generation of `.chem` files given a molecule's 3D structure (e.g. PDB file).
        - Its goal is to enable the automatization of SMIFs pipelines for ligands.
        - Currently only implemented for stacking interactions, following a heuristic geometrical approach.
    - `pwoverlap` subcommand for calculation of point-wise overlaps between two biomolecular structures.
        - Motivation: SMIF grid approach can be sensitive to resolution/coordinate system of the grid when evaluating SMIF values at precise points in space.
        - A **point-wise (PW)** overlap is proposed instead, where the same underlying model from SMIFs is used, but at concrete points in space instead of a discrete grid.


## [0.11.1] - 2026-05-12
- Added option to run smiffer on both protein and nucleic atoms at the same time, by using the `macro` subcommand (in place of `prot` or `rna`).
- Fixed bug in `smutils occupancy`.
- Fixed bug where only one stacking group was considered per residue name from the chemtable files.
- Set default *Occupancy Grid* radii back to 2.0.


## [0.11.0] - 2026-05-07
- Important refactorings in some internal implementations (mostly for smiffer and the APBS wrapper).
    - Improved memory usage.
    - Moved `log_apbs` to an independent smutils tool.
    - Removed *SMIF_HYDRODIFF*, as it is now unnecesary.
        - This operation can now be performed very easily as postprocessing by using `vgtools op sub` (e.g. `volgrids vgtools op sub hphob.mrc hphil.mrc hdiff.mrc`)
    - All SMIF and trimming operations are now ensured to use the same molecular system information.
        - This can be either the input structure directly, or use atom data from the pertinent intermediary PQR file (in case `DO_SMIF_APBS=true`).
    - Other internal refactorings.

- Added `vgtools points` for getting the grid values at specific points in space.
- Trajectory grids now base their box on the global minimum/maximum coordinates accross all frames (instead of just the first one).
- Standard DNA residue names (**DT**, **DC**, **DA**, **DG**) and alternative RNA residue names (**RU**, **RC**, **RA**, **RG**) should now be recognized by smiffer.
    - **DT** will be considered as analogous to **U** for SMIF operations.


## [0.10.0] - 2026-05-04
- General:
    - Trajectory mode can now be handled with optional multiprocessing (`-n` flag).
    - `GridIO.read_auto` now accepts an optional "key" argument (to be used when reading CMAP files).
    - `GridIO.write_auto` now takes an optional key argument for saving CMAP.
    - OG radius config can now be individually set for every OG kind.

- Bug fixes:
    - Fixed bug in `ResiduesNucleic` (bad parsing of rnapolis results).
    - Fixed bug where `GridIO.get_cmap_keys()` was crashing if its input file didn't exist.
    - `_warning_big_grid()` was incorrectly placed as a method for `Box`, when it should be a method for `Grid`.

- Vendors:
    - Updated vendor `freyacli` to v0.4.0.
    - Updated vendor `molutils` to v0.4.0.
    - **NOTE**: These updates break compatibility with previous versions.
    - Added automatic generation of a `versions.txt` file to register the intended vendor versions for a given volgrids version.

- VGTools:
    - Internal `average` and `rotate` methods now return the grid instead of saving it.
        - The saving is performed inside `AppVGTools`. That way, the operation methods can be used by other scripts without the intermediary saving step.
    - Split internal `op` method into two: unary and binary operations.
        - Both internal methods are now iterators that yield the resulting grids (CMAP case) instead of saving them.


## [0.9.0] - 2026-04-29
- General:
    - Fixed compatibility issues with older versions of Python and rnapolis.
    - Changed default configuration to `OVERWRITE_OK=true` until evaluating if the overwrite safeguard is properly implemented (mostly in terms of convenience).
    - Vendor packages now consider only the local copy.
        - It should be clear now that they're a light dependency and don't need to be pip installed.

- New features `vgtools`:
    - Added `AND` and `OR` operations to `vgtools op`.
        - This operation chooses the points between two grids according to which one has the maximum absolute value.
    - Added option to `vgtools op` (with the flag `-b`) for interpolating to smallest enclosing box.
        - This is useful for for keeping all the points between the two grids involved in an operation, in case they aren't perfectly overlaping (e.g. one completely enclosing the other one).

- Residue calculations (`smiffer`, `smutils`):
    - Residue-specific calculations now follow the syntax **chain_id.resid** (insted of simply **resid**).
    - Changed flags and subcommands to reflect this new behavior:
        - `smiffer`'s `-i/--resids` changed to `-r/--residues`.
        - `smutils`'s `resids_nonbp` changed to `res_nobp`.
    - Path to a list with residues is no longer accepted as a value for the `-r` flag.
    - Added `res_nostk` to `smutils`.
        - `res_nostk` is similar to `res_nobp`: it selects residues that are not already participating in at least 2 stacking interactions.
        - This will also work for DNA structures (thanks to a temporary fix in the RNA .chem tables).

- For future use:
    - Added new `KOperation` `ABSMAX` for internal usage (eventually).
    - Added AMBER charges and radius data, as used by `pdb2pqr`
        - This might be useful later on for a possible way of implementing `OgAPBS`, or with any other application that needs those partial charges.


## [0.8.2] - 2026-04-24
- Occupancy grids:
    - OGs are now stamped with OR operation (instead of addition)
    - Stacking OG now considers a sphere for every stacking atom, instead of a single sphere in the center of geometry
- Operations from `vgtools op` can now handle CMAP files with multiple grids.
- Attempting to overwrite existing file now asks for confirmation.
    - This behavior can be disabled with by setting the config `OVERWRITE_OK=true`
- Added `smutils histogram` utility.
- Updated packaged freyacli vendor to version `0.3.0`.


## [0.8.1] - 2026-04-22
- Started implementing new features relevant for an *interfaces* pipeline.
    - Added a generic `op` to `vgtools`.
        - **op** stands for arithmetic *operation*. Currently implemented are: `abs`, `add`, `sub`, `mul`, `div`
    - Added `occupancy` grids to `smutil`.
- Reorganized the files containing the SMIF classes.


## [0.8.0] - 2026-04-21
- Massive refactoring of the CLI logic.
    - Delegating the CLI parsing and management to the [freyacli](https://github.com/DiegoBarMor/freyacli) package.
        - Help strings are no longer awkwardly hardcoded, but automatically generated by freyacli instead.
        - Argument rules should be more robust now too.
    - Added some colors to the displayed outputs via the `freyacli.Color` methods.
    - Simplified the CLI for some of the `vgtools`.
- Freyacli is included as a vendor when distributing to pip, so it's not necessary to install it.
    - If volgrids isn't installed via pip, then freyacli is automatically fetched the first time volgrids is executed.


## [0.7.0] - 2026-04-16
- New features:
    - Started the **Smiffer Utilities** application. Use it via `volgrids smutils`.
        - First utility: `volgrids smutils resids_nonpb`.
            - Use it to print a space-separated list of RNA residue indices for nucleotides that are not currently engaged in a canonical base-pairing interaction.
            - Requires **rnapolis** to be installed (`pip install rnapolis`).
    - Added `--resids` (`-i`) flag to **Smiffer**.
        - A custom list of residue indices can now be specified. SMIFs will be calculated for only those residues.
        - Trimming and APBS remain unaffected.

- Tweaks:
    - Renamed `DO_SIMPLE_HBONDS_RNA` to `HBONDS_ONLY_NUCLEOBASE`.
    - An error is now raised for invalid configuration keys.
    - A list of known configuration keys can be printed by passing the -c flag with no values following.
        - It still needs to have the required arguments behind. So use it like e.g. `volgrids smiffer rna 1akx.pdb -c`.


## [0.6.1] - 2026-04-10
- Cavities can now also be used to *weight* the smifs instead of trimming them.
    - The `CAVITIES_WEIGHT` configuration controls how much the smifs inside cavities should be taken into account.
    - Positive `CAVITIES_WEIGHT` values give more weight to smifs inside cavities; negative values give more weight to smifs outside.
- Added `DO_SIMPLE_HBONDS_RNA` configuration (default: false)
    - When set to true, the RNA hbonds will only consider atoms from the nucleobases.
    - It can remove a significant amount of noise from the grid, at the cost of some relevant interactions lost.


## [0.6.0] - 2026-04-10
- Added a `vgtools` utility for rotating grids.
- Implemented config option for ensuring *equilateral grids*.
    - These are grids that have the same amount of points in every dimension i.e. an homogeneous resolution. This is independent from the deltas.
    - Equilateral grids are naturally heavier, but they could be convenient for some operations.
- Improved the naive cavities finder by adding the possibility to perform multiple passes in evenly spaced rotations (controlled by the `CAVITIES_NPASSES` configuration)
- The `--config` (`-c`) flag can now also take a list of either config files or *keyword=value* pairs, or a mix of both.
- the **apbs** application and the `apbs.sh` utility now displayer errors properly.
- Some general code improvements.


## [0.5.1] - 2026-04-07
- INI parser now considers the standard `;` as comment start character.
    - This affects the `config.ini` files for the parameters and the `.chem` tables.


## [0.5.0] - 2026-04-05
- Implemented experimental **cavities** trimming operation.
- Added a `vgtools` utility for summarizing grids.


## [0.4.2] - 2026-04-03
- Bugfixes:
    - Fixed trimming mask being saved with an incorrect name.
    - When installing with pip, required packages should now be installed automatically corectly.
- New QOL features:
    - Volgrids version is now displayed with e.g. `volgrids --help` (or simply `volgrids`).
    - Old cmap output is now cleared when running smiffer.


## [0.4.1] - 2026-01-27
- Fixed bug where APBS grids weren't trimmed properly.
    - Cause: in the `APBS` pipeline, an intermediary `PQR` file is generated by `pdb2pqr`, which can add extra atoms in an attempt to fix the structures. Such atoms can include whole functional groups that weren't being considered before by `smiffer`'s trimmers.
- The intermediary `PQR` file is now saved as an output when dealing with a single structure (i.e. outside trajectory mode).


## [0.4.0] - 2026-01-16
- Bug fixes:
    - Fixed `config_volgrids.ini` never being used.
    - Fixed `ParserIni` parsing commented headers as uncommented.
- `ParserIni` now captures the text preceding the first header into an empty string key.
- Configuration files no longer use headers (headers can still be placed in the files for clarity, but they will be inconsequential).
- Added a separate timer for APBS subprocesses.


## [0.3.0] - 2026-01-15
- Fixed compatibility issue with Python 3.14.
- Added `GridIO.write_auto` method for generic output operations (format infered from path's suffix).
- Package requirements are now automatically checked/installed when installing VolGrids with PIP.
- Changed the way VolGrids should be run from the root directory of VolGrid's repo (when not installed). Instead of having separate entry scripts for every "application" (e.g. `smiffer.py`, `vgtools.py`...), these can be accessed by running `python3 volgrids` followed by the application name as an argument. For example, `python3 volgrids smiffer`.
- Similarly, if VolGrids is installed as a package, it can be called as an independent command from anywere. The previous example would correspond to `volgrids smiffer` instead.
- Added an **apbs** app to VolGrids. It intends to simplify the process of running APBS on structures without having to worry about intermediary files. Run by `volgrids apbs`.
- Changed the behavior of the `-a` flag for `smiffer`. Once again, it's only valid if an APBS output files follows; in such a case, than this file is used to generate the electrostatic SMIF. If the flag is not used, the application falls back to automatically calculating the APBS output as a temporary file.


## [0.2.0] - 2026-01-13
- APBS can now be automatically computed by Smiffer.
- Electrostatic smifs (APBS) can now be also generated in trajectory mode.


## [0.1.0] - 2026-01-13
- Initial upload to PyPI.
