import numpy as np
from pathlib import Path

import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class NumberLists:
    def __init__(self, data: list[str], width: int) -> None:
        """
        `data` is provided as a flat list of strings.
        It can be either a list of space-separated numbers, or the path to CSV file with comma-separated numbers [WIP].
        `width` is the number of values per row of the parsed numbers table. It also corresponds to the columns expected in the CSV file.
        If there is an error in the parsing, `NumberLists.error` will contain the error message, and `NumberLists.values` will be empty.
        Otherwise, `NumberLists.values` will contain the parsed numbers as a tuple of tuples, where each inner tuple corresponds to a row of the parsed numbers table.
        """
        self.error: str | None = None
        self.values: tuple[tuple[float|int]] | np.ndarray = ()
        self._parse(data, width)


    # --------------------------------------------------------------------------
    def _parse(self, data: list[str], width: int) -> None:
        if len(data) == 1 and Path(data[0]).is_file():
            self._from_csv(data[0], width)
            return

        if len(data) % width:
            self.error = f"Got {len(data)} values, which is not a multiple of {width}."
            return

        self.values = tuple(zip(*(
            map(float, data[i::width])
            for i in range(width)
        )))


    # --------------------------------------------------------------------------
    def _from_csv(self, path_csv: Path, width: int) -> None:
        """
        Read a CSV of with numerical values. It must have `width` columns.
        An optional non-numeric header row is allowed and silently ignored.
        """
        rows = np.atleast_2d(np.genfromtxt(path_csv, delimiter = ",", dtype = vg.FLOAT_DTYPE))

        ##### genfromtxt parses a non-numeric header row as NaN; drop any such rows
        self.values = rows[~np.isnan(rows).any(axis = 1)]
        if self.values.size == 0: raise ValueError(
            f"No valid numeric rows found in box CSV file: {path_csv}"
        )
        if self.values.shape[1] != 6: raise ValueError(
            f"CSV file must have exactly {width} columns "
            f"but got {self.values.shape[1]} columns in: {path_csv}"
        )


# //////////////////////////////////////////////////////////////////////////////
