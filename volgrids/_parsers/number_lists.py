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
        self.values: tuple[tuple[float|int]] = ()
        self._parse(data, width)


    # --------------------------------------------------------------------------
    def _parse(self, data: list[str], width: int) -> None:
        if len(data) % width:
            self.error = f"Got {len(data)} values, which is not a multiple of {width}."
            return

        self.values = tuple(zip(*(data[i::width] for i in range(width))))


# //////////////////////////////////////////////////////////////////////////////
