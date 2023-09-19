"""
Unit tests for open reading frame.
"""
from sequence_analysis.open_reading_frame import OpenReadingFrame

class TestOpenReadingFrame:
    def test_init(self):
        orf = OpenReadingFrame()
        orf.start = 5
        assert orf.start == 5
        orf.stop = 55
        assert orf.stop == 55
        orf.parent_sequence = "ACGT"
        assert orf.parent_sequence == "ACGT"