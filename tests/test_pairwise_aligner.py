"""
Unit tests for pairwise aligner.
"""
from sequence_analysis import PairwiseAligner

import pdb

class TestPairwiseAligner:
    def test_init(self):
        pa = PairwiseAligner()
        pa.query = "ACGT"
        pa.target = "ACCT"
        
        assert pa.algorithm == "global"
        assert pa.query == "ACGT"
        assert pa.target == "ACCT"

    def test_init_2(self):
        pa = PairwiseAligner(query="ACGT", target="ACCT", algorithm="global")
        assert pa.algorithm == "global"
        assert pa.query == "ACGT"
        assert pa.target == "ACCT"

    def test_align(self):
        pa = PairwiseAligner()
        pa.query = "HEAGAWGHEE"
        pa.target = "PAWHEAE"

        pa.align()

    def test_align_2(self):
        pa = PairwiseAligner()
        pa.algorithm = "local"
        pa.query = "HEAGAWGHEE"
        pa.target = "PAWHEAE"

        pa.set_gap_penalty(8)

        pa.align()

    def test_align_blastn(self):
        pa = PairwiseAligner("blastn")
        pa.algorithm = "local"

        pa.target = "GGTACGTACG"
        pa.query = "ACCT"

        pa.align()