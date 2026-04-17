# TODO

<!-- ----------------------------------------------------------------------- -->
## General
* improve help strings. Add annotations, docstrings and overall cleaning
* replace mdanalysis with another PDB parser?
    * find out how to handle trajectory files in this case
* improve the CLI parsing
    * e.g. passing another positional arguments after `volgrids smiffer rna path_input.pdb` shouldn't be allowed


<!-- ----------------------------------------------------------------------- -->
## VolGrids
* implement: raise an error if a format file is opened with the wrong function
* add tests for parameters being directly passed to the App classes (instead of parsing the CLI arguments)


<!-- ----------------------------------------------------------------------- -->
## SMIFFER
* list of spheres for trajectory+pocket_sphere mode
* check what happens with structure files with multiple models.
* maybe: replace the RNDS trimming with a faster method
* change the ligand example to one that uses both NAMES_HBACCEPTORS, NAMES_HBDONORS and NAMES_HBD_FIXED
* document the .chem tables
* check if there's a bug in the peptide bond N of the test toy system peptide_no_h
* add safeguard when there's no atoms for the specified molecule type
* add tests for apbs
* add tests for -i flag
* add possibility for treshold i.e. removing low value points (treshold of 0.5 already can reduce CMAP sizes by 90%)
* reimplement automatic script generation for visualizing pockets in VMD (pocket-sphere mode)
* check whether the trimmer can still be saved when no smifs are calculated
* reduce the number of file writes performed during trajectory mode
  * introduce multithreading for calculating multiple frames concurrently?


<!-- ----------------------------------------------------------------------- -->
## VEINS
* finish/rework "energies" mode implementation
* implement "forces" mode
* move Grid's static fields into config_volgrids.ini
* add tests


<!-- ----------------------------------------------------------------------- -->
## VGTools
* check what happens if performing "fix_cmap" operation when cmap input and output are the same file
* implement the fixing operation directy on "packing", to ensure that packed frames have the same resolution (add flag to override this behavior)
* mode to perform operations on grids: abs, sum, diff, mask...
* when editing a CMAP file (be it converting it or performing an operation on it), one should be able to specify the key of the relevant grid (instead of GridIO.read_auto arbitrarily deciding to take the first key it finds in the CMAP header)
* bypass the "large grid" warning when processing an existing large grid with VGTools.
* add tests for the "average", "summary", "rotate" operations.
* give a warning if "convert" doesn't take any output option.
* improve CMAP keys automatically assigned when packing.


<!-- ----------------------------------------------------------------------- -->
## SMUtils
* example operations:
    * `findsphere` to find the position and (optionally extendable) radius of a sphere surrounding a query for the structure.
    * `smifhist` for an interactive matplotlib visualization of the SMIFs' histograms.


<!-- ----------------------------------------------------------------------- -->
