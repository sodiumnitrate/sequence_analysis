"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files, and various filtering methods.
"""
import sequence_analysis_cpp

class SeqSet(sequence_analysis_cpp.SeqSet):
    """This class holds a list of sequence objects of a given type."""
    def __str__(self):
        """__str__ function for sequence set (seq_set) object."""
        return f"Sequence set object with {len(self)} sequences of type {self.type}"