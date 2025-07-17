import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class ParserConfig(vg.ParserIni):
    def apply_config(self,
        section: str, scope_module: dict[str, ], scope_dependencies: dict[str, ],
        valid_config_keys: set[str], all_configs_mandatory: bool = True
    ) -> None:
        """
        Applies the configuration to the provided global dictionary.
        """

        if not self.has(section):
            raise ValueError(f"Configuration file does not contain [{section}] section (case sensitive).")

        for k, value in self.iter_splitted_lines(section):
            k = k.upper()
            if k not in valid_config_keys:
                raise ValueError(f"Unknown configuration for [{section}]: {k}.")
            scope_module[k.upper()] = self._parse_str(scope_dependencies, value)
            valid_config_keys.remove(k)

        if all_configs_mandatory and valid_config_keys:
            raise ValueError(f"Configuration keys not set: {', '.join(valid_config_keys)}.")


    # --------------------------------------------------------------------------
    @staticmethod
    def _parse_str(scope: dict, str_value: str):
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
        stripped = str_value.strip('"')
        if len(stripped) != len(str_value):
            return stripped

        stripped = str_value.strip("'")
        if len(stripped) != len(str_value):
            return stripped

        ### MODULE ATTRIBUTES
        error_message = f"Invalid configuration value: {str_value}"

        parts = str_value.split('.')
        if len(parts) < 2: raise ValueError(error_message)

        key = parts.pop(0)
        node = scope.get(key, None)
        if node is None: raise ValueError(error_message)

        while parts:
            key = parts.pop(0)
            try:
                node = node.__dict__.get(key, None)
            except (AttributeError):
                raise ValueError(error_message)
            if node is None:
                raise ValueError(error_message)
        return node


# //////////////////////////////////////////////////////////////////////////////
