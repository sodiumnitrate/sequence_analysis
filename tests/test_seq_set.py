"""
Basic tests for the seq_seq class.
"""
from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence

class TestSeqSet:
    def test_set_type(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')
        assert sequences.type == 'dna'

        sequences = seq_set(file_name="aux_files/pub.fasta")
        assert sequences.type == 'protein'

    def test_slice_1(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')

        seq = sequences[0]

        assert isinstance(seq, sequence)

        assert seq == sequences.records[0]

    def test_slice_2(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')

        new_set = sequences[:2]

        assert isinstance(new_set, seq_set)
        assert len(new_set) == 2