import re

# //////////////////////////////////////////////////////////////////////////////
class ParserIni:
    """Parser for INI-like files, extracting sections and their corresponding lines."""
    _CHAR_COMMENT = '#'

    # --------------------------------------------------------------------------
    def __init__(self, path_ini):
        self._ini_sections: dict[str, list[str]] = {}
        with open(path_ini, 'r') as file:
            self._extract_sections(file.read())


    # --------------------------------------------------------------------------
    @classmethod
    def split_line(cls, line: str, sep: str = '=') -> tuple[str, str]:
        """
        Splits a line into a key-value pair based on the specified separator.
        """
        line = line.split(cls._CHAR_COMMENT)[0].strip() # Remove comments
        pair = tuple(map(str.strip, line.split(sep)))
        if len(pair) != 2:
            raise ValueError(f"Line '{line}' does not contain '{sep}'")
        return pair


    # --------------------------------------------------------------------------
    def has(self, key: str) -> bool:
        """
        Checks if the given key exists in the sections.
        """
        return key in self._ini_sections.keys()


    # --------------------------------------------------------------------------
    def get(self, key, default = None) -> list[str]:
        """
        Returns the value for the given key, or default if the key is not found.
        """
        return self._ini_sections.get(key, default)


    # --------------------------------------------------------------------------
    def iter_lines(self, key: str):
        """
        Iterates over the lines of the section identified by key.
        Yields each line that is not empty or a comment.
        """
        for line in self._ini_sections.get(key, []):
            line = line.strip()
            if not line or line.startswith(self._CHAR_COMMENT): continue
            yield line


    # --------------------------------------------------------------------------
    def iter_splitted_lines(self, key: str, sep: str = '='):
        """
        Iterates over the lines of the section identified by key,
        splitting each line into a key-value pair based on the specified separator.
        Yields tuples of (key, value) for each line.
        """
        for line in self.iter_lines(key):
            yield self.split_line(line, sep)


    # --------------------------------------------------------------------------
    def _extract_sections(self, data: str) -> None:
        """
        Extracts all headers and their corresponding bodies from the given string.
        Each header is a substring enclosed in square brackets, e.g. [HEADER].
        The body is the text after the header's closing bracket up to the next header or end of string.
        """
        pattern = re.compile(r"\[(\w+)\]")
        matches = list(pattern.finditer(data))
        for i, match in enumerate(matches):
            header = match.group(1)
            body_start = match.end()
            body_end = matches[i + 1].start() if (i + 1 < len(matches)) else len(data)
            body = data[body_start:body_end].strip()
            if header in self._ini_sections: raise ValueError(f"Duplicate header found: {header}")
            self._ini_sections[header] = body.splitlines()


# //////////////////////////////////////////////////////////////////////////////
