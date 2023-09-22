"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files, and various filtering methods.
"""
from .sequence_analysis_cpp import SeqSet as SeqSet_cpp

class SeqSet(SeqSet_cpp):
    """This class holds a list of sequence objects of a given type."""
    def __str__(self):
        """__str__ function for sequence set (seq_set) object."""
        return f"Sequence set object with {len(self)} sequences of type {self.type}"

    # TODO: overload constructor to enable initialization with a list of sequences and automatic type setting.
