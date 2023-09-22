"""
Unit tests for GenomeMap.
"""
from sequence_analysis import GenomeMap

class TestGenomeMap:
    def test_init(self):
        gm = GenomeMap()
        gm.chromosome_name = "CH2"
        assert gm.chromosome_name == "CH2"

        gm.sample_name = "tissue"
        assert gm.sample_name == "tissue"
