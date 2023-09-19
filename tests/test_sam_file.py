"""
Unit tests for the SamFile class.
"""
from sequence_analysis.sam_file import SamFile

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