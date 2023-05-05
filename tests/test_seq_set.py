"""
Basic tests for the seq_seq class.
"""
from sequence_analysis.seq_set import seq_set

class TestSeqSet:
    def test_set_type(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')
        assert sequences.type == 'dna'

        sequences = seq_set(file_name="aux_files/pub.fasta")
        assert sequences.type == 'protein'