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

    def test_init_3(self):
        pa = PairwiseAligner("blastn")

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

    def test_align_3(self):
        pa = PairwiseAligner()
        pa.algorithm = "local"
        pa.query = "MPA"
        pa.target = "MKA"
        pa.align()

    def test_get_alignment_range(self):
        pa = PairwiseAligner()
        pa.algorithm = "local"
        pa.query = "MCDDVAALVVD"
        pa.target = "DVAALV"
        pa.align()

        qs, qe = pa.get_query_range()
        assert qs == 3
        assert qe == 8

        ts, te = pa.get_target_range()
        assert ts == 0
        assert te == len(pa.target)-1

    def test_align_blastn(self):
        pa = PairwiseAligner("blastn", -2)
        pa.algorithm = "local"

        pa.target = "GGTACGTACG"
        pa.query = "ACCT"

        pa.align()

    def test_levenshtein(self):
        pa = PairwiseAligner("levenshtein")
        pa.algorithm = "local"

        pa.query = "AAACGTGG"
        pa.target = "AACGTT"
        pa.align()

        assert len(pa.get_query_aligned()) == len(pa.get_target_aligned())
        assert len(pa.get_match_string()) == len(pa.get_query_aligned())
        assert pa.get_score() == 5

        assert '.' not in pa.get_match_string()

    def test_levenshtein_2(self):
        pa = PairwiseAligner("levenshtein")
        pa.algorithm = "global"

        pa.query = "AAACGTGG"
        pa.target = "AACGTT"
        pa.align()

        assert len(pa.get_query_aligned()) == len(pa.get_target_aligned())
        assert len(pa.get_match_string()) == len(pa.get_query_aligned())
        assert pa.get_score() == 5

        assert '.' in pa.get_match_string()

        assert pa.get_alignment_start() == 0
        assert pa.get_alignment_end() == len(pa.target)

    def test_levenshtein_3(self):
        pa = PairwiseAligner("levenshtein")
        pa.algorithm = "local"

        pa.query = "ACT"
        pa.target = "CGACTGAC"
        pa.align()

        assert pa.get_alignment_start() == 2
        assert pa.get_alignment_end() == 5
