"""
Unit tests for StarMap

TODO: how to properly test for this?
"""
from sequence_analysis.star_map import StarMap

class TestStarMap:
    def test_init(self):
        sm = StarMap(out_folder="aux_files")
        assert sm.n_threads == 1