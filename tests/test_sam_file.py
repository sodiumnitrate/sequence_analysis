"""
Unit tests for the SamFile class.
"""
import numpy as np

from sequence_analysis import SamFile

class TestSamFile:
    def test_init(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        assert sf.file_name is not None
        assert sf.get_normalized()

    def test_set_filter(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0], [-1], ["CM042140.1"], 0)

    def test_read_1(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0], [-1], ["CM042140.1"], 0)
        sf.read()
        assert len(sf) == 10

        sf2 = SamFile()
        sf2.file_name = "aux_files/sam_test.sam"
        sf2.set_filter_options([0], [-1], ["CM042177.1"], 0)
        sf2.read()
        assert len(sf2) == 2

    def test_read_2(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read()
        assert sf.get_normalized()
        assert len(sf) == 9

        assert sf.get_seq_start() != 0
        assert sf.get_seq_end() == 128727393

    def test_read_3(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.read()
        assert len(sf) == 12

    def test_read_4(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0,0], [-1,-1], ["CM042140.1", "CM042177.1"], 0)
        sf.read()
        assert len(sf) == 12

    def test_normalize(self):
        sf = SamFile()
        sf.file_name = "aux_files/reference_out.sam"
        sf.set_normalized_false()
        assert not sf.get_normalized()
        sf.get_lengths_from_fasta("aux_files/reference_set_dna.fasta")
        sf.read()

    def test_get_genome_map(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read()

        gm = sf.get_genome_map("CM042140.1", "tissue")

        assert gm.sample_name == "tissue"
        assert gm.chromosome_name == "CM042140.1"

        heatmap = gm.get_heatmap(sf.get_seq_start(), sf.get_seq_end())

        assert isinstance(heatmap, list)
        assert np.sum(heatmap) > 0

        heatmap = gm.get_heatmap(sf.get_seq_start()+100, sf.get_seq_start()+200)

        assert len(heatmap) == 101

    def test_sam_file_add(self):
        sam1 = SamFile()
        sam1.file_name = "aux_files/sam_test.sam"
        sam1.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sam1.read()

        sam2 = SamFile()
        sam2.file_name = "aux_files/sam_test.sam"
        sam2.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sam2.read()

        sam1.add_sam_file(sam2)
        assert len(sam1) == 2*len(sam2)

    def test_get_entries_headers(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0], [-1], ["CM042140.1"], 97)
        sf.read()

        assert len(sf) == 9
        entries = sf.get_entries()
        assert len(entries) == 9
        for el in entries:
            assert "CM042140.1" in el

        headers = sf.get_headers()
        assert len(headers) == 1
        assert "CM042140.1" in headers[0]

    def test_get_entries_headers_2(self):
        sf = SamFile()
        sf.file_name = "aux_files/sam_test.sam"
        sf.set_filter_options([0,0], [-1,-1], ["CM042140.1", "CM042177.1"], 0)
        sf.read()
        assert len(sf) == 12
        entries = sf.get_entries()
        for el in entries:
            assert "CM042140.1" in el or "CM042177.1" in el

        headers = sf.get_headers()
        assert len(headers) == 2
        for el in headers:
            assert "CM042140.1" in el or "CM042177.1" in el
