"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files, and various filtering methods.
"""
from pathlib import Path
from .sequence_analysis_cpp import SeqSet as SeqSet_cpp

class SeqSet(SeqSet_cpp):
    """This class holds a list of sequence objects of a given type."""
    def __init__(self, file_name=None, list_of_seqs=None):
        """
        Overload the constructor so that you can set file name in one go.

        TODO: do this on the C++ side so that you don't have to worry about keyword arguments?
        """
        # init object
        super(SeqSet, self).__init__()

        # check if file_name or list_of_seqs provided
        if file_name is not None:
            if list_of_seqs is not None:
                print("ERROR: you provided both a file_name and list_of_seqs. Not setting anything.")
                raise ValueError
            if Path(file_name).suffix in [".fasta", ".fna", ".fa", ".ffn", ".frn", ".fst"]:
                self.read_fasta(file_name)
            else:
                print(f"ERROR: file with suffix {Path(file_name).suffix} not supported.")
                raise ValueError
        
        if list_of_seqs is not None:
            if not isinstance(list_of_seqs, list):
                raise TypeError

            self.records = list_of_seqs

    def __str__(self):
        """__str__ function for sequence set (seq_set) object."""
        return f"Sequence set object with {len(self)} sequences of type {self.type}"

    def add_records(self, list_of_sequences):
        if not isinstance(list_of_sequences, list):
            raise TypeError

        for el in list_of_sequences:
            self.add_sequence(el)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.records[index]
        elif isinstance(index, slice):
            new_set = SeqSet()
            new_set.add_records(self.records[index])
            return new_set