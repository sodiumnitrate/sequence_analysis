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
    def __init__(self, seq_str, name=None):
        """
        Overload the constructor on C++ side so that you can auto-assign type.
        """
        super(Sequence, self).__init__(seq_str)
        self.set_type()
        if name is not None:
            if not isinstance(name, str):
                raise TypeError
            self.name = name

    def __getitem__(self, key):
        """Overwrites __getitem__ so that a part of the sequence can be returned."""
        if isinstance(key, slice) or isinstance(key, int):
            return self.seq_str[key]
        else:
            raise TypeError

    def select_dna_as_protein(self, protein_sequence):
        """
        If we have a dna sequence, select a subsequence that corresponds
        to the given protein sequence.

        TODO: make seq object a valid input
        """
        if not isinstance(protein_sequence, str):
            raise TypeError

        if self.type != 'dna' and self.type != 'rna':
            raise TypeError

        orfs = self.get_open_reading_frames()
        for orf in orfs:
            idx = orf.protein_sequence.find(protein_sequence)
            if idx != -1:
                seq = Sequence(orf.rna_sequence[idx*3:idx*3+len(protein_sequence)*3])
                return seq.seq_str

        return None
