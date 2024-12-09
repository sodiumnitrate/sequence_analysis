"""
Unit tests for the SamFile class.
"""
import numpy as np
import sys
sys.path.insert(0, '.') 

from sequence_analysis import SamFile

class TestSamFile:
    def test_init(self):
        sf = SamFile()
        assert sf.get_normalized()

    def test_set_filter(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 0)

    def test_read_1(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 0)
        sf.read("aux_files/sam_test.sam")
        assert len(sf) == 10

        sf2 = SamFile()
        sf2.set_filter_options([0], [-1], ["CM042177.1"], 0)
        sf2.read("aux_files/sam_test.sam")
        assert len(sf2) == 2

    def test_read_2(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read("aux_files/sam_test.sam")
        assert sf.get_normalized()
        assert len(sf) == 9

        assert sf.get_seq_start() != 0
        assert sf.get_seq_end() == 128727393

    def test_read_3(self):
        sf = SamFile()
        sf.read("aux_files/sam_test.sam")
        assert len(sf) == 12

    def test_read_4(self):
        sf = SamFile()
        sf.set_filter_options([0,0], [-1,-1], ["CM042140.1", "CM042177.1"], 0)
        sf.read("aux_files/sam_test.sam")
        assert len(sf) == 12

    def test_normalize(self):
        sf = SamFile()
        sf.set_normalized_false()
        assert not sf.get_normalized()
        sf.get_lengths_from_fasta("aux_files/reference_set_dna.fasta")
        sf.read("aux_files/reference_out.sam")

        scores = sf.get_normalized_scores()
        for score in scores:
            assert score >= 0 and score <= 100

    def test_get_genome_map(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read("aux_files/sam_test.sam")

        gm = sf.get_genome_map("CM042140.1", "tissue")

        assert gm.sample_name == "tissue"
        assert gm.chromosome_name == "CM042140.1"

        heatmap = gm.get_heatmap(sf.get_seq_start(), sf.get_seq_end())

        assert isinstance(heatmap, list)
        assert np.sum(heatmap) > 0

        heatmap = gm.get_heatmap(sf.get_seq_start()+100, sf.get_seq_start()+200)

        assert len(heatmap) == 101

        mapped_names = gm.get_mapped_read_names(sf.get_seq_start(), sf.get_seq_end())
        assert len(mapped_names) == 9

    def test_get_genome_map_multi(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read("aux_files/sam_test.sam")

        # double counts
        multiplicity = {"SRR18071790.134779":2,
                        "SRR18071790.134809":2,
                        "SRR18071790.134866":2,
                        "SRR18071790.134932":2,
                        "SRR18071790.135003":2,
                        "SRR18071790.135067":2,
                        "SRR18071790.135115":2,
                        "SRR18071790.135306":2,
                        "SRR18071790.135320":2,
                        "SRR18071790.135420":6}

        sf.set_multiplicity(multiplicity)
        gm = sf.get_genome_map("CM042140.1", "tissue")
        heatmap1 = gm.get_heatmap(sf.get_seq_start()+100, sf.get_seq_start()+200)

        sf1 = SamFile()
        sf1.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf1.read("aux_files/sam_test.sam")

        gm1 = sf1.get_genome_map("CM042140.1", "tissue")
        heatmap2 = gm1.get_heatmap(sf.get_seq_start()+100, sf.get_seq_start()+200)

        assert np.sum(heatmap1) == 2 * np.sum(heatmap2)


    def test_sam_file_add(self):
        sam1 = SamFile()
        sam1.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sam1.read("aux_files/sam_test.sam")

        sam2 = SamFile()
        sam2.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sam2.read("aux_files/sam_test.sam")

        sam1.add_sam_file(sam2)
        assert len(sam1) == 2*len(sam2)

    def test_sam_file_add_2(self):
        sam1 = SamFile()
        sam1.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sam1.read("aux_files/sam_test.sam")

        sam2 = SamFile()
        sam2.copy_filters_from_another(sam1)
        sam2.add_sam_file(sam1)

        assert len(sam1) == len(sam2)

    def test_get_headers(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read("aux_files/sam_test.sam")

        assert len(sf) == 9

        headers = sf.get_headers()
        assert len(headers) == 1
        assert "CM042140.1" in headers[0]

    def test_get_headers_2(self):
        sf = SamFile()
        sf.set_filter_options([0,0], [-1,-1], ["CM042140.1", "CM042177.1"], 0)
        sf.read("aux_files/sam_test.sam")
        assert len(sf) == 12

        headers = sf.get_headers()
        assert len(headers) == 2
        for el in headers:
            assert "CM042140.1" in el or "CM042177.1" in el

    def test_copy_filters(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        starts = sf.get_starts()
        ends = sf.get_ends()
        names = sf.get_names()
        sf.read("aux_files/sam_test.sam")

        sf2 = SamFile()
        sf2.copy_filters_from_another(sf)
        starts2 = sf.get_starts()
        ends2 = sf.get_ends()
        names2 = sf.get_names()

        assert starts == starts2
        assert ends == ends2
        assert names == names2

        sf2.read("aux_files/sam_test.sam")

        assert len(sf) == len(sf2)

    def test_read_parallel(self):
        sf = SamFile()
        sf.read_multiple_files(["aux_files/sam_test.sam", "aux_files/sam_test_2.sam"])

        assert len(sf) == 24

    def test_get_multimapping_stats(self):
        sf = SamFile()
        sf.read("aux_files/sam_test.sam")

        multimappers_dict = sf.get_multimapping_stats()
        assert len(multimappers_dict['SRR18071790.135420']) == 3

    def test_to_seq_set(self):
        sf = SamFile()
        sf.read("aux_files/sam_test.sam")

        sset = sf.to_seq_set()
        assert len(sset) == len(sf)
        sset.remove_duplicates()
        assert len(sset) == len(sf) - 2
        assert sset.records[0].name.startswith("SRR180")

    def test_get_names(self):
        sf = SamFile()
        sf.read("aux_files/sam_test.sam")

        names = sf.get_seq_names()
        assert isinstance(names, set)
        assert len(names) == 10

    def test_read_bin_flag_filter(self):
        sf = SamFile()
        sf.set_primary_map_true()
        sf.read("aux_files/sam_test.sam")
        assert len(sf.get_entries()) == 11

    def test_read_bin_flag_filter_2(self):
        sf = SamFile()
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.set_primary_map_true()
        sf.read("aux_files/sam_test.sam")
        assert len(sf.get_entries()) == 8