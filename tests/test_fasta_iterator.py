"""
Unit tests for the fasta iterator object.
"""
from sequence_analysis import FastaIterator
from sequence_analysis import Sequence
import pdb

class TestFastaIterator:
    def test_init(self):
        fasta = FastaIterator("aux_files/new_test.fasta")

    def test_read(self):
        fasta = FastaIterator("aux_files/new_test.fasta")
        seq = fasta.get_next()

        assert seq is not None
        assert isinstance(seq.name, str)
        assert seq.name == ' seq0'
        assert isinstance(seq.seq_str, str)

    def test_read_2(self):
        fasta = FastaIterator("aux_files/new_test.fasta")
        ct = 1
        seq = fasta.read()
        while seq:
            seq = fasta.read()
            ct += 1

        assert ct == 4