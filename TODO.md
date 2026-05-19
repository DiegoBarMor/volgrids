# TODO

<!-- ----------------------------------------------------------------------- -->
## General
* Add annotations, docstrings and overall cleaning
* Start replacing mdanalysis with another PDB parser
    * its overhead can be avoided if just simple PDB parsing is required.
        * keep using MDAnalysis for trajectory files, as well as fallback for other molecular formats.
    * remove the dependency to `gridData` in `grid_io.py`
        * use the [mrcfile](https://mrcfile.readthedocs.io/en/stable/) library directly for handling MRC/CCP4 files.
        * DX parser (read) could be manually implemented.
* improve the calls to `AppMain.load_configs`?


<!-- ----------------------------------------------------------------------- -->
## VolGrids
* add explanations to the list of configs printed with the empty `-c` flag
* implement: raise an error if a format file is opened with the wrong function
* add tests for parameters being directly passed to the App classes (instead of parsing the CLI arguments)
* check if the implementation of the OVERWRITE_OK flag is user-convenient


<!-- ----------------------------------------------------------------------- -->
## SMIFFER
* remove the need to write a prot/rna/macro/ligand subcommand.
* can't use the -r flag when PQR structure is used because it doesn't have chain information.
* check what happens with structure files with multiple models.
* list of spheres for trajectory+pocket_sphere mode
* maybe: replace the RNDS trimming with a faster method
* change the ligand example to one that uses both NAMES_HBACCEPTORS, NAMES_HBDONORS and NAMES_HBD_FIXED
* document the .chem tables
* check if there's a bug in the peptide bond N of the test toy system peptide_no_h
* add safeguard when there's no atoms for the specified molecule type
* add tests for apbs
* add tests for -r flag
* add possibility for treshold i.e. removing low value points (treshold of 0.5 already can reduce CMAP sizes by 90%)
* reimplement automatic script generation for visualizing pockets in VMD (pocket-sphere mode)
* check whether the trimmer can still be saved when no smifs are calculated


<!-- ----------------------------------------------------------------------- -->
## VEINS
* finish/rework "energies" mode implementation
* implement "forces" mode
* move Grid's static fields into config_volgrids.ini
* add tests


<!-- ----------------------------------------------------------------------- -->
## VGTools
* add missing vgtools tests.
* check what happens if performing "fix_cmap" operation when cmap input and output are the same file
* implement the fixing operation directy on "packing", to ensure that packed frames have the same resolution (add flag to override this behavior)
* when editing a CMAP file (be it converting it or performing an operation on it), one should be able to specify the key of the relevant grid (instead of GridIO.read_auto arbitrarily deciding to take the first key it finds in the CMAP header)
* bypass the "large grid" warning when processing an existing large grid with VGTools.
* give a warning if "convert" doesn't take any output option.
* improve the CMAP keys that are automatically assigned when packing.


<!-- ----------------------------------------------------------------------- -->
## SMUtils
* check why `DEBUG_CHEMTABLE_LIGAND` config isn't being taken into account
* ideas:
    * `findsphere` to find the position and (optionally extendable) radius of a sphere surrounding a query for the structure.


<!-- ----------------------------------------------------------------------- -->
