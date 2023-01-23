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

    def test_choose_in_between_matching(self):
        seq1 = sequence("AAAAMDMAAAAAMDMAAA")
        matching, spans, inverted = seq1.choose_all_matching_patterns("MDM",return_between_matching=True)
        assert(seq1.seq[spans[0][0]:spans[0][1]+1] == "MDM")
        assert(seq1.seq[inverted[0][0]:inverted[0][1]+1] == "AAAAA")

    def test_six_frame(self):
        seq = sequence("gcgctgaaagcgctgattatggatatggcgctgaaagcgctgatt".upper())
        check = seq.six_frame_check("MDM")
        # TODO: why is there an X at the end?
        assert(check == "ALKALIMDMALKALI")

        check2 = seq.six_frame_check("WWWW")
        assert(check2 == None)