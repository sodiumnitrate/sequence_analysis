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

    def test_align(self):
        pa = PairwiseAligner()
        pa.query = "HEAGAWGHEE"
        pa.target = "PAWHEAE"

        pa.align()