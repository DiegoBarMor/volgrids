import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class ConfigParser(vg.IniParser):
    def apply_config(self,
        key: str, scope: dict[str, ], valid_configs = set[str],
        all_configs_mandatory: bool = True
    ) -> None:
        """
        Applies the configuration to the provided global dictionary.
        """

        if not self.has(key):
            raise ValueError(f"Configuration file does not contain [{key}] section (case sensitive).")

        for k, value in self.iter_splitted_lines(key):
            k = k.upper()
            if k not in valid_configs:
                raise ValueError(f"Unknown configuration for [{key}]: {k}.")
            scope[k.upper()] = self._parse_str(scope, value)
            valid_configs.remove(k)

        if all_configs_mandatory and valid_configs:
            raise ValueError(f"Configuration keys not set: {', '.join(valid_configs)}.")


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
