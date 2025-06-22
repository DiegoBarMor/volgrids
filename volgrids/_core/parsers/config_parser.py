import volgrids as vg

# //////////////////////////////////////////////////////////////////////////////
class ConfigParser(vg.IniParser):
    def apply_config(self,
        key: str, globs: dict[str, ], valid_configs = set[str],
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
            globs[k.upper()] = self.parse_str(value)
            valid_configs.remove(k)

        if all_configs_mandatory and valid_configs:
            raise ValueError(f"Configuration keys not set: {', '.join(valid_configs)}.")


# //////////////////////////////////////////////////////////////////////////////
