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

    def test_init_from_list(self):
        gm = GenomeMap()
        gm.set_from_list([0, 5, 15, 20], [10, 8, 25, 22])
        heatmap = gm.get_heatmap(0,25)

        assert len(heatmap) == 26
        assert heatmap[2] == 1
        assert heatmap[5] == 2
        assert heatmap[8] == 2
        assert heatmap[9] == 1
        assert heatmap[12] == 0
        assert heatmap[25] == 1