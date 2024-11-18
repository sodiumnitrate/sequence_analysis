"""
Python bindings for the SamEntry object.
"""
from .sequence_analysis_cpp import SamEntry as SamEntry_cpp

class SamEntry(SamEntry_cpp):
    """
    Python bindings for the SamEntry class.
    """
    pass
