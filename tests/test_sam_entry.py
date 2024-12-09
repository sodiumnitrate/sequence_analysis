"""
Unit tests for the SamEntry class.
"""
from sequence_analysis import SamFile, SamEntry
import sys
sys.path.insert(0, '.')

class TestSamEntry:
    def test_init(self):
        entry = SamEntry('b', 'b', 'b', 1, 2,  0.1, 0)

    def test_read_and_extract(self):
        sf = SamFile()
        sf.set_filter_options([0,0], [-1,-1], ["CM042140.1", "CM042177.1"], 0)
        sf.read("aux_files/sam_test.sam")

        entries = sf.get_entries()

        assert len(entries) == len(sf)

    def test_get_props(self):
        sf = SamFile()
        sf.set_filter_options([0,0], [-1,-1], ["CM042140.1", "CM042177.1"], 0)
        sf.read("aux_files/sam_test.sam")

        entries = sf.get_entries()

        seq = entries[0].get_seq_str()
        rn = entries[0].get_read_name()
        m = entries[0].get_mapped_onto()
        
        assert entries[0].get_length() == len(entries[0])

    def test_get_binary_flag_0(self):
        entry = SamEntry('b', 'b', 'b', 1, 2,  0.1, 256)
        assert entry.get_binary_flag() == 256

    def test_get_binary_flag_1(self):
        sf = SamFile()
        sf.read('aux_files/sam_single.sam')
        entry = sf.get_entries()[0]
        assert entry.get_binary_flag() == 256

    def test_get_binary_flag_2(self):
        sf = SamFile()
        sf.read('aux_files/sam_test.sam')
        entries = sf.get_entries()

        assert len(entries) == 12

        assert entries[0].get_binary_flag() == 0
        assert entries[6].get_binary_flag() == 256
