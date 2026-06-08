# Backend

Any operations that are written in a language different then Python (or bash) are to be placed in subdirectories here.
They will be called (and compiled, if necessary) by the pertinent Python scripts.

- Current backend code:
    - `segmentation`:
        - Given a volgrid's **input grid**, segment the regions of space enclosed by a certain **isovalue** (as visualized in VMD, ChimeraX, etc).
        - Clusters under a certain **volume threshold** (that is, a minimum accepted number of points inside the cluster) are discarded.
        - Clusters are sorted in descending order according to their volume. They are labelled with integer values in the **output grid**.
        - Both the input and output grids must be in `BIN` format.
