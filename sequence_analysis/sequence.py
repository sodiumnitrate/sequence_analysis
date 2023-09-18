"""
This file holds the sequence class and related methods.
"""
import sequence_analysis_cpp

class Sequence(sequence_analysis_cpp.Sequence):
    """
    Python bindings for the Sequence class. Inherits the C++ class to add
    a bunch of python-specific functions.
    """
    def __getitem__(self, key):
        """Overwrites __getitem__ so that a part of the sequence can be returned."""
        if isinstance(key, slice) or isinstance(key, int):
            return self.seq_str[key]
        else:
            raise TypeError