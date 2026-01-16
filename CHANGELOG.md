# Changelog

## [0.4.0] - 2026-01-16
- Bug fixes:
    - Fixed `config_volgrids.ini` never being used.
    - Fixed `ParserIni` parsing commented headers as uncommented.
- `ParserIni` now captures the text preceding the first header into an empty string key.
- Configuration files no longer use headers (headers can still be placed in the files for clarity, but they will be inconsequential).
- Added a separate timer for APBS subprocesses.

## [0.3.0] - 2026-01-15
- Fixed compatibility issue with Python 3.14.
- Added `GridIO.write_auto` method for generic output operations (format infered from path's suffix).
- Package requirements are now automatically checked/installed when installing VolGrids with PIP.
- Changed the way VolGrids should be run from the root directory of VolGrid's repo (when not installed). Instead of having separate entry scripts for every "application" (e.g. `smiffer.py`, `vgtools.py`...), these can be accessed by running `python3 volgrids` followed by the application name as an argument. For example, `python3 volgrids smiffer`.
- Similarly, if VolGrids is installed as a package, it can be called as an independent command from anywere. The previous example would correspond to `volgrids smiffer` instead.
- Added an **apbs** app to VolGrids. It intends to simplify the process of running APBS on structures without having to worry about intermediary files. Run by `volgrids apbs`.
- Changed the behavior of the `-a` flag for `smiffer`. Once again, it's only valid if an APBS output files follows; in such a case, than this file is used to generate the electrostatic SMIF. If the flag is not used, the application falls back to automatically calculating the APBS output as a temporary file.

## [0.2.0] - 2026-01-13
- APBS can now be automatically computed by Smiffer.
- Electrostatic smifs (APBS) can now be also generated in trajectory mode.

## [0.1.0] - 2026-01-13
- Initial upload to PyPI.
