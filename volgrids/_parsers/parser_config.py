import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class ParserConfig(vg.ParserIni):
    def apply_config(self,
        scope_module: dict[str,],
        this_module_keys: set[str],
    ) -> None:
        """
        Applies the configuration to the provided global dictionary.
        """
        keys = tuple(self._ini_sections.keys())
        if len(keys) > 1:
            raise ValueError(
                f"Configuration file must have no headers, but the following were found: {keys[1:]}. "+\
                "Remove all headers and keep only key-value pairs (or comments) in the file."
            )
        assert keys[0] == '', f"Unexpected INI key with value '{keys[0]}'. Check ParserIni."

        for k, value in self.iter_splitted_lines(''):
            k = k.upper()
            if k not in vg.KNOWN_CONFIGS: raise ValueError(f"Unknown configuration: {k}.")
            if k not in this_module_keys: continue
            scope_module[k] = self._parse_str(value)


    # --------------------------------------------------------------------------
    @staticmethod
    def _parse_str(str_value: str):
        ### INTEGERS
        if str_value.isdigit():
            return int(str_value)

        ### FLOATS
        try: return float(str_value)
        except ValueError: pass

        ### BOOLEANS
        if str_value.lower() in ["true", "false"]:
            return str_value.lower() == "true"

        ### STRINGS
        return str_value.strip('"').strip("'")


# //////////////////////////////////////////////////////////////////////////////
