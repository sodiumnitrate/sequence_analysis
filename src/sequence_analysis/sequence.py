"""
This file holds the sequence class and related methods.
"""
from .sequence_analysis_cpp import Sequence as Sequence_cpp

class Sequence(Sequence_cpp):
    """
    Python bindings for the Sequence class. Inherits the C++ class to add
    a bunch of python-specific functions.

    std::string name;
    std::string seq_str;
    std::string type;
    """
    def __init__(self, seq_str):
        """
        Overload the constructor on C++ side so that you can auto-assign type.
        """
        super(Sequence, self).__init__(seq_str)
        self.set_type()

    def __getitem__(self, key):
        """Overwrites __getitem__ so that a part of the sequence can be returned."""
        if isinstance(key, slice) or isinstance(key, int):
            return self.seq_str[key]
        else:
            raise TypeError
