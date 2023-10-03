"""
Unit tests for StarMap

TODO: how to properly test for this?
"""
import sys
sys.path.insert(0, '.') 

from sequence_analysis import StarMap

class TestStarMap:
    def test_init(self):
        sm = StarMap(out_folder="aux_files")
        assert sm.n_threads == 1
