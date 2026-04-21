### [TODO] the overlap methods will most likely be reconsidered as generic "op" (add, sub, mul, div, abs) subcommands

# # --------------------------------------------------------------------------
# def _parse_overlap(self) -> None:
#     self._set_help_str(
#         "usage: volgrids vgtools overlap [grid1] [grid2] [output] [options...]",
#         "Compute overlap between two molecular interaction fields.",
#         "The first grid will be interpolated to match the coordinate system of the second grid.",
#         "Available options:",
#         "    -h, --help         Show this help message and exit.",
#         "    --op, --operation  Operation type: 'multiply' (default), 'subtract', 'add'.",
#     )
#     if self._has_param_kwds("help"):
#         self._exit_with_help()

#     vgt.PATH_OVERLAP_GRID1 = self._safe_path_file_in(
#         self._safe_get_param_pos(1,
#            err_msg = "No first grid file provided. Provide a path to the first grid file as second positional argument."
#         )
#     )

#     vgt.PATH_OVERLAP_GRID2 = self._safe_path_file_in(
#         self._safe_get_param_pos(2,
#            err_msg = "No second grid file provided. Provide a path to the second grid file as third positional argument."
#         )
#     )

#     vgt.PATH_OVERLAP_OUT = self._safe_path_file_out(
#         self._safe_get_param_pos(3,
#            err_msg = "No output grid file provided. Provide a path where to save the overlap grid as fourth positional argument."
#         )
#     )

#     vgt.OVERLAP_OPERATION = self._safe_kwd_str("operation", default="multiply")


# # --------------------------------------------------------------------------
# def _parse_overlap_cross(self) -> None:
#     vgt.PATH_OVERLAP_CROSS_GRID1 = self._safe_path_file_in(
#         self._safe_get_param_pos(1,
#            err_msg = "No first grid file provided. Provide a path to the first grid file as second positional argument."
#         )
#     )

#     vgt.PATH_OVERLAP_CROSS_GRID2 = self._safe_path_file_in(
#         self._safe_get_param_pos(2,
#            err_msg = "No second grid file provided. Provide a path to the second grid file as third positional argument."
#         )
#     )

#     vgt.PATH_OVERLAP_CROSS_OUT = self._safe_path_file_out(
#         self._safe_get_param_pos(3,
#            err_msg = "No output grid file provided. Provide a path where to save the overlap grid as fourth positional argument."
#         )
#     )


# # --------------------------------------------------------------------------
# def _parse_overlap_diff(self) -> None:
#     vgt.PATH_OVERLAP_DIFF_GRID1 = self._safe_path_file_in(
#         self._safe_get_param_pos(1,
#            err_msg = "No first grid file provided. Provide a path to the first grid file as second positional argument."
#         )
#     )

#     vgt.PATH_OVERLAP_DIFF_GRID2 = self._safe_path_file_in(
#         self._safe_get_param_pos(2,
#            err_msg = "No second grid file provided. Provide a path to the second grid file as third positional argument."
#         )
#     )

#     vgt.PATH_OVERLAP_DIFF_OUT = self._safe_path_file_out(
#         self._safe_get_param_pos(3,
#            err_msg = "No output grid file provided. Provide a path where to save the difference grid as fourth positional argument."
#         )
#     )


# //////////////////////////////////////////////////////////////////////////////
