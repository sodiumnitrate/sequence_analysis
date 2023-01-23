from hashlib import new
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

    def test_needleman_wunsch_3(self):
        seq1 = "HEAGAWGHEE"
        seq2 = "PAWHEAE"

        new_alignment = pairwise_alignment(seq1, seq2)
        new_alignment.align()

        x = new_alignment.sequence1_aligned.seq
        y = new_alignment.sequence2_aligned.seq

        assert(new_alignment.score == 1)
        assert(x == "HEAGAWGHE-E")
        assert(y == "--P-AW-HEAE")

    def test_smith_waterman_1(self):
        seq1 = "HEAGAWGHEE"
        seq2 = "PAWHEAE"
        new_alignment = pairwise_alignment(seq1,seq2,algorithm="smith-waterman")
        new_alignment.align()

        x = new_alignment.sequence1_aligned.seq
        y = new_alignment.sequence2_aligned.seq

        assert(new_alignment.score == 28)
        assert(x == "AWGHE")
        assert(y == "AW-HE")

    def test_smith_waterman_2(self):
        seq1 = "ALKALIMDMALKALI"
        seq2 = "WWWWWWMDMKALI"
        new_alignment = pairwise_alignment(seq1,seq2,algorithm="smith-waterman")
        new_alignment.align()

        x = new_alignment.sequence1_aligned.seq
        y = new_alignment.sequence2_aligned.seq

        assert(new_alignment.score == 29)
        assert(x == "MALKALI")
        assert(y == "MDMKALI")

    # TODO: check empty string or string of numbers, etc.
    # pytest coverage?