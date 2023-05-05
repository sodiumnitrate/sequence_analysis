from matplotlib import use
from sequence_analysis.pairwise_alignment import pairwise_alignment

class TestPairwiseAlignment:
    def test_needleman_wunsch_1(self):
        seq1 = "ALKALI"
        seq2 = "ALKALI"
        new_alignment = pairwise_alignment(seq1, seq2)
        new_alignment.use_blosum_50 = True
        new_alignment.align()
        x = new_alignment.target_aligned.seq
        y = new_alignment.query_aligned.seq
        assert (x == y)

    def test_needleman_wunsch_2(self):
        seq1 = "ALKALIMDMKALI"
        seq2 = "ALKAKIWMKAL"
        new_alignment = pairwise_alignment(seq1, seq2)
        new_alignment.use_blosum_50 = True
        new_alignment.align()

        x = new_alignment.target_aligned.seq
        y = new_alignment.query_aligned.seq
        assert (new_alignment.score == 29)
        assert (x == "ALKALIMDMKALI")
        assert (y == "ALKAKIW-MKAL-")

    def test_needleman_wunsch_3(self):
        seq1 = "HEAGAWGHEE"
        seq2 = "PAWHEAE"

        new_alignment = pairwise_alignment(seq1, seq2)
        new_alignment.use_blosum_50 = True
        new_alignment.align()

        x = new_alignment.target_aligned.seq
        y = new_alignment.query_aligned.seq

        assert (new_alignment.score == 1)
        assert (x == "HEAGAWGHE-E")
        assert (y == "--P-AW-HEAE")

    def test_smith_waterman_1(self):
        seq1 = "HEAGAWGHEE"
        seq2 = "PAWHEAE"
        new_alignment = pairwise_alignment(
            seq1, seq2, algorithm="smith-waterman")
        new_alignment.use_blosum_50 = True
        new_alignment.align()

        x = new_alignment.target_aligned.seq
        y = new_alignment.query_aligned.seq

        assert (new_alignment.score == 28)
        assert (x == "AWGHE")
        assert (y == "AW-HE")

    def test_smith_waterman_2(self):
        seq1 = "ALKALIMDMALKALI"
        seq2 = "WWWWWWMDMKALI"
        new_alignment = pairwise_alignment(
            seq1, seq2, algorithm="smith-waterman")
        new_alignment.use_blosum_50 = True
        new_alignment.align()

        x = new_alignment.target_aligned.seq
        y = new_alignment.query_aligned.seq

        assert (new_alignment.score == 29)
        assert (x == "MALKALI")
        assert (y == "MDMKALI")

    def test_biopython_vs_needleman_wunsch(self):
        seq1 = "HEAGAWGHEE"
        seq2 = "PAWHEAE"
        alignment_nw = pairwise_alignment(
            seq1, seq2, algorithm="needleman-wunsch", gap=0, gap_open=0)
        alignment_nw.use_blosum_50 = True
        alignment_nw.align()
        alignment_bp = pairwise_alignment(
            seq1, seq2, algorithm="biopython-global")
        alignment_bp.use_blosum_50 = True
        alignment_bp.align()

        assert (alignment_nw.score == alignment_bp.score)
        # NOTE: in general, it's unreasonable to expect the following to match because of the ambiguity in situations where two scores are equal while forming F_ij. One possibility is to have the option to return all possible alignments, from biopython's alignment.
        # assert(alignment_nw.sequence1_aligned.seq == alignment_bp.sequence1_aligned.seq)

    def test_biopython_vs_needleman_wunsch_custom_params(self):
        seq1 = "HEAGAWGHEE"
        seq2 = "PAWHEAE"
        alignment_nw = pairwise_alignment(
            seq1,
            seq2,
            match=2,
            unmatch=-1,
            gap=-0.1,
            gap_open=-0.1,
            algorithm="needleman-wunsch")
        alignment_nw.align()
        alignment_bp = pairwise_alignment(
            seq1,
            seq2,
            match=2,
            unmatch=-1,
            gap=-0.1,
            gap_open=-0.1,
            algorithm="biopython-global")
        alignment_bp.align()

        assert (alignment_nw.score == alignment_bp.score)
        # TODO: extend function to access all possibilities within the list of alignments biopython returns?
        # assert(alignment_nw.sequence1_aligned.seq == alignment_bp.sequence1_aligned.seq)

    # TODO: check empty string or string of numbers, etc.

    def test_biopython_local(self):
        # TODO: check robustness w.r.t. biopython version
        seq1 = "TACCG"
        seq2 = "ACG"

        alignment = pairwise_alignment(
            seq1,
            seq2,
            match=1,
            unmatch=0,
            gap=0,
            gap_open=0,
            algorithm="biopython-local")
        alignment.align()

        assert alignment.target_aligned.seq == "ACCG"
        assert alignment.query_aligned.seq == "AC-G"
