from sequence_analysis.pairwise_alignment import pairwise_alignment

class TestPairwiseAlignment:
    def test_needleman_wunsch_1(self):
        seq1 = "ALKALI"
        seq2 = "ALKALI"
        new_alignment = pairwise_alignment(seq1,seq2)
        new_alignment.align()
        x = new_alignment.sequence1_aligned.seq
        y = new_alignment.sequence2_aligned.seq
        assert(x == y)

    def test_needleman_wunsch_2(self):
        seq1 = "ALKALIMDMKALI"
        seq2 = "ALKAKIWMKAL"
        new_alignment = pairwise_alignment(seq1,seq2)
        new_alignment.align()

        x = new_alignment.sequence1_aligned.seq
        y = new_alignment.sequence2_aligned.seq
        assert(new_alignment.score == 29)
        assert(x == "ALKALIMDMKALI")
        assert(y == "ALKAKIW-MKAL-")