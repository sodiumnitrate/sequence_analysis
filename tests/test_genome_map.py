"""
Unit tests for GenomeMap.
"""
from sequence_analysis import GenomeMap

import pdb

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

    def test_init_from_list_2(self):
        gm = GenomeMap()
        gm.set_from_list([0, 5, 15, 20], [10, 8, 25, 22], [2, 1, 1, 2])
        heatmap = gm.get_heatmap(0,25)

        assert len(heatmap) == 26
        assert heatmap[2] == 2
        assert heatmap[5] == 3
        assert heatmap[8] == 3
        assert heatmap[9] == 2
        assert heatmap[12] == 0
        assert heatmap[20] == 3
        assert heatmap[25] == 1

    def test_init_from_list_3(self):
        gm = GenomeMap()
        gm.set_from_list([0, 5, 15, 20], [10, 8, 25, 22], [0, 0, 0, 0])
        heatmap = gm.get_heatmap(0, 25)

        assert len(heatmap) == 26
        assert all([val == 0 for val in heatmap])

    def test_init_from_list_4(self):
        gm = GenomeMap()
        starts = [112523435, 111365804, 110177291, 109994207, 109207913, 107708375, 106422323, 105030476, 109376135, 110133209, 111060220, 110354704, 110116726, 109964725, 107209270, 105222459, 107192592, 107572674, 109980828, 110193438, 110255241, 110028084, 109721307, 108676134, 110014822]

        ends = [112525600, 111367945, 110180182, 109996783, 109210078, 107710537, 106424491, 105032602, 109378246, 110136082, 111062310, 110356809, 110119437, 109967586, 107211387, 105224549, 107194898, 107574854, 109983746, 110196290, 110257364, 110030984, 109723466, 108678407, 110017740]

        mult = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        gm.set_from_list(starts, ends, mult)

        gm2 = GenomeMap()
        gm2.set_from_list(starts, ends)