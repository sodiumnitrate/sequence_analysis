from sequence_analysis.sequence import sequence

class TestSequence:
    def test_choose_matching(self):
        seq1 = sequence("AAAAMDMAAAAAMDMAAA")
        matching, spans = seq1.choose_all_matching_patterns("MDM")
        assert(matching[0] == "MDM")
        assert(len(matching)==2)

        assert(spans[0][0] == 4)
        assert(spans[0][1] == 6)
        assert(spans[1][0] == 12)
        assert(spans[1][1] == 14)
