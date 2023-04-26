from sequence_analysis.sequence import sequence
import pdb

class TestIndexWithoutGaps:
    def test_with_gaps(self):
        seq = sequence('---AAA---TTT')
        index = seq.find_index_with_gaps(index_with_gaps=3)
        assert index == 0

    def test_without_gaps(self):
        seq = sequence('---A-AA-TT-T')
        index = seq.find_index_with_gaps(index_without_gaps=0)
        assert index == 3

    def test_with_gaps_2(self):
        seq = sequence('AA-AA-TTT-AA--')
        index = seq.find_index_with_gaps(index_with_gaps=2)
        assert index is None

        index = seq.find_index_with_gaps()
        assert index is None

        index = seq.find_index_with_gaps(index_with_gaps=1)
        assert index == 1

        index = seq.find_index_with_gaps(index_with_gaps=4)
        assert index == 3

        index = seq.find_index_with_gaps(index_without_gaps=5)
        assert index == 7