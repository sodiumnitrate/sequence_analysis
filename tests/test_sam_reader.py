"""
Unit tests for sam_reader.
"""
from sequence_analysis.sam_reader import SamReader
import pdb


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

    def test_normalize(self):
        sr = SamReader("aux_files/out.sam")
        sr.read()

        sr.read_full_lengths("aux_files/ubiquitin.fasta")
        assert len(sr.sam_string_list) == 2

        assert len(sr.full_lengths) == 1

        sr.normalize_AS_and_filter(min_score=50)
        assert len(sr.sam_string_list) == 0
        
        sr.read()
        sr.normalize_AS_and_filter(min_score=1)
        assert len(sr.sam_string_list) == 2

    def test_get_unique(self):
        sr = SamReader("aux_files/sam_test.sam")
        sr.read()

        ref_names = sr.get_unique_reference_names()
        assert len(ref_names) == 2
        