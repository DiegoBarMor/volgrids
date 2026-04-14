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
        for section in self._ini_sections.keys():
            for k, value in self.iter_splitted_lines(section):
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
