"""
Unit tests for genome map.
"""
import numpy as np
from sequence_analysis.genome_map import GenomeMap
import pdb

class TestGenomeMap:
    def test_init(self):
        gm = GenomeMap('CHNAME', 3)

        assert len(gm.file_names) == 0
        assert gm.ch_name == 'CHNAME'
        assert gm.ch_num == 3
        assert gm.sample_name is None

        assert gm.start == 0
        assert gm.end == np.inf

        assert gm.min_AS == 0

        assert len(gm.reads) == 0
        assert gm.read_range == (0, np.inf)

    def test_filters(self):
        gm = GenomeMap('CHNAME', 3)
        gm.set_alignment_score_filter(88)
        gm.set_range_filter(10,5000)

        assert gm.start == 10
        assert gm.end == 5000
        assert gm.read_range == (10, 5000)
        assert gm.min_AS == 88

    def test_get_reads(self):
        gm = GenomeMap('CM042140.1', 4)
        gm.set_alignment_score_filter(98)
        gm.file_names = ['aux_files/sam_test.sam']
        gm.get_reads()

        assert len(gm.reads) > 0

    def test_get_heatmap(self):
        gm = GenomeMap('CM042140.1', 4)
        gm.set_alignment_score_filter(98)
        gm.file_names = ['aux_files/sam_test.sam']

        gm.get_reads()
        print(gm.read_range)

        gm.get_heatmap()
        assert len(gm.heatmap) == gm.read_range[1] - gm.read_range[0]
        assert np.sum(gm.heatmap.freqs) > 0

    def test_get_heatmap_2(self):
        gm = GenomeMap('CM042140.1', 4)
        gm.set_alignment_score_filter(98)
        gm.file_names = ['aux_files/sam_test.sam'] * 2
        gm.get_reads()

        gm.get_heatmap()
        assert len(gm.heatmap) == gm.read_range[1] - gm.read_range[0]
        assert np.sum(gm.heatmap.freqs) > 0