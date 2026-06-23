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
* refactor some of the inheritance in SMIF classes (specially for HBonds).
* improve the calls to `AppMain.load_configs`?


<!-- ----------------------------------------------------------------------- -->
## VolGrids
* generalize the usage of the `-c` flag (for customazing configurations) in all modes.
* add explanations to the list of configs printed with the empty `-c` flag
* implement: raise an error if a format file is opened with the wrong function
* add tests for parameters being directly passed to the App classes (instead of parsing the CLI arguments)
* check if the implementation of the OVERWRITE_OK flag is user-convenient
* idea: centralize the usage of `mda.Universe` instances into a single wrapper class (MolSystem is already there) and add `delete_traj_locks` in its destructor (will it work?).
* deal somehow with using a wrong comment char inside a config file


<!-- ----------------------------------------------------------------------- -->
## SMIFFER
* validate the INI headers when parsing a user's provided chem table
* there seems to be some previous bugs in trajectory mode; track down and fix.
* check what happens with structure files with multiple models.
* RNDS trimming could be removed. instead of it, a post-processing with a similar (and probably better) effect could be done with SMIF segmentation.
* change the ligand example to one that uses both HBACCEPTORS, HBDONORS and NAMES_HBD_FIXED
* document the .chem tables
* check if there's a bug in the peptide bond N of the test toy system peptide_no_h
* add tests for apbs
* add tests for -r flag
* add possibility for treshold i.e. removing low value points (treshold of 0.5 already can reduce CMAP sizes by 90%)
* check whether the trimmer can still be saved when no smifs are calculated.
* MDAnalysis struggles when parsing certain PQR files.


<!-- ----------------------------------------------------------------------- -->
## VGTools
* add missing vgtools tests.
* check what happens if performing "fix_cmap" operation when cmap input and output are the same file
* implement the fixing operation directy on "packing", to ensure that packed frames have the same resolution (add flag to override this behavior)
* when editing a CMAP file (be it converting it or performing an operation on it), one should be able to specify the key of the relevant grid (instead of `Grid.load` arbitrarily deciding to take the first key it finds in the CMAP header)
* bypass the "large grid" warning when processing an existing large grid with VGTools.
* add tests for BIN format and for "segment".
* let "vgtools summary" receive multiple input files.


<!-- ----------------------------------------------------------------------- -->
## SMUtils
* check why `DEBUG_CHEMTABLE_LIGAND` config isn't being taken into account
* reimplement automatic script generation for visualizing pockets in VMD (pocket-sphere mode)


<!-- ----------------------------------------------------------------------- -->
