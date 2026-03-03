from __future__ import annotations

from ._core import __doc__, __version__,parse_frame_from_file

__all__ = ["__doc__", "__version__","parse_frame_from_file"]

# # __init__.py
# # import the compiled extension module (name discovered from wheel)
# from . import lammpsParser as _ext  # if compiled module is lammpsParser.so
# # or if compiled is _lammpsparser:
# # from . import _lammpsparser as _ext

# # re-export conveniences
# parse_frame_from_file = _ext.parse_frame_from_file
# parse_frame_to_numpy = _ext.parse_frame_to_numpy
# parse_frames_from_directory = _ext.parse_frames_from_directory
# install_into_active_env = _ext.install_into_active_env

# # expose raw extension if needed
# _extmod = _ext
