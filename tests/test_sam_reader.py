"""
Unit tests for sam_reader.
"""
from sequence_analysis.sam_reader import SamReader

class TestSamReader:
    def test_sam_reader_1(self):
        sr = SamReader("aux_files/sam_test.sam", mapped_onto="CM042177.1")
        sr.read()

        assert len(sr.sam_string_list) == 2
        assert len(sr.header) == 1

    def test_sam_reader_2(self):
        sr = SamReader("aux_files/sam_test.sam")
        sr.read()

        assert len(sr.sam_string_list) == 12
        assert len(sr.header) == 2

    def test_sam_reader_3(self):
        sr = SamReader("aux_files/sam_test.sam", start=115000000)
        sr.read()

        assert len(sr.sam_string_list) == 7

    def test_sam_reader_4(self):
        sr = SamReader("aux_files/out.sam")
        sr.read()
        assert len(sr.sam_string_list) == 2