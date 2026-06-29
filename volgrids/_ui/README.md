# UI

Volgrids uses [FreyaCLI](https://github.com/DiegoBarMor/freyacli) for managing its CLI. This directory contains very specific files related the UI:
- `app_main.py`: Entry point for the volgrids suite of applications.
- `app_subcommand.py`: Base class for all subcommands of the volgrids suite of applications.
- `config_manager.py`: Class for handling of configs used by volgrids and its sub-applications.
- `fy_rules.fyr`: `FYR` file that freyacli uses to determine the allowed CLI inputs received from a user.
- `fy_help.fyh`: `FYH` file that freyacli uses to display help strings when requested / parsing a bad input.
