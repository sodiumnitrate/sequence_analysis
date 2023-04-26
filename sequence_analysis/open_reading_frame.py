"""
This file holds the open reading frame class.
"""
import sequence_analysis.sequence as seq

class OpenReadingFrame:
    """
    This class holds the open reading frame information for a given dna or 
    rna sequence.
    """

    def __init__(self,
                 rna_sequence,
                 parent_sequence,
                 start,
                 stop,
                 strand,
                 frame,
                 protein_sequence=None):

        self.rna_sequence = rna_sequence
        self.parent_sequence = parent_sequence
        self.start = start
        self.stop = stop
        self.strand = strand
        self.frame = frame
        self.protein_sequence = protein_sequence

        self.check_init_data()
        if protein_sequence is None:
            self.protein_sequence = self.rna_sequence.translate()

    def check_init_data(self):
        assert isinstance(self.rna_sequence, seq.sequence)
        assert isinstance(self.parent_sequence, seq.sequence)
        assert isinstance(self.start, int)
        assert isinstance(self.stop, int)
        assert isinstance(self.frame, int)
        assert self.frame == 0 or self.frame == 1 or self.frame == 2
        assert isinstance(self.protein_sequence, seq.sequence) or self.protein_sequence is None

    
