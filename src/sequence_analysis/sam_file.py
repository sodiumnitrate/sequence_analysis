"""
Python bindings for the SamFile object.
"""
from .sequence_analysis_cpp import SamFile as SamFile_cpp

class SamFile(SamFile_cpp):
    """
    Python bindings for the SamFile class.
    """
    pass
