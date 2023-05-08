from sequence_analysis.sequence import sequence
from sequence_analysis.seq_set import seq_set

import pdb

class TestSequence:
    def test_translate(self):
        seq = sequence("ATGGATCGGTTCGGTAGGGCTCGATCACATCGCTAG")
        assert seq.type == 'dna'
        seq = seq.translate()
        assert seq.type == 'protein'
        assert seq.seq == 'MDRFGRARSHR*'

    def test_translate_2(self):
        seq = sequence("ATGGATC--GGTTC----GGTAGGGCTCGATCACATCGCTAG")
        assert seq.type == 'dna'
        seq = seq.translate()
        assert seq.type == 'protein'
        assert seq.seq == 'MDRFGRARSHR*'

    def test_reverse_complement_rna(self):
        seq = sequence("AUGAUGC")
        assert seq.type == 'rna'
        seq2 = seq.reverse_complement()
        assert seq2.seq == "GCAUCAU"

    def test_choose_matching(self):
        seq1 = sequence("AAAAMDMAAAAAMDMAAA")
        matching, spans = seq1.choose_all_matching_patterns("MDM")
        assert (matching[0] == "MDM")
        assert (len(matching) == 2)

        assert (spans[0][0] == 4)
        assert (spans[0][1] == 6)
        assert (spans[1][0] == 12)
        assert (spans[1][1] == 14)

    def test_choose_in_between_matching(self):
        seq1 = sequence("AAAAMDMAAAAAMDMAAA")
        _, spans, inverted = seq1.choose_all_matching_patterns(
            "MDM", return_between_matching=True)
        assert (seq1.seq[spans[0][0]:spans[0][1] + 1] == "MDM")
        assert (seq1.seq[inverted[0][0]:inverted[0][1] + 1] == "AAAAA")

    def test_frame_shift(self):
        seq = sequence("gcgctgaaagcgctgattatggatatggcgctgaaagcgctgatt".upper())
        seq2 = seq.frame_shift(frame=0)
        assert (len(seq2) % 3 == 0)
        seq2 = seq.frame_shift(frame=1)
        assert (len(seq2) % 3 == 0)
        seq2 = seq.frame_shift(frame=2)
        assert (len(seq2) % 3 == 0)

    def test_frame_shift_2(self):
        seq = sequence("gcgctgaaagcgctgattatggatatggcgctgaaagcgctgatt".upper())
        seq2 = seq.frame_shift(frame=0)
        assert len(seq2) == len(seq)
        seq2 = seq.frame_shift(frame=1)
        assert len(seq2) == len(seq)
        seq2 = seq.frame_shift(frame=2)
        assert len(seq2) == len(seq)

    def test_six_frame(self):
        seq = sequence("atggcgctgaaagcgctgattatggatatggcgctgaaagcgctgatttag".upper())
        check, _ = seq.six_frame_check("MDM", min_orf_len=0)
        assert (check == "MALKALIMDMALKALI")

        check2 = seq.six_frame_check("WWWW")
        assert (check2[0] is None)
        assert (check2[1] is None)

    def test_choose_matching_2(self):
        set1 = seq_set(file_name="aux_files/pub.fasta")
        a1 = set1.records[0]
        matching, spans, inverted = a1.choose_all_matching_patterns(
            "(?:PER[WYQH])*(?:[MF]D(?:(?![MF]D).){4,6}){3,10}", return_between_matching=True)

        assert (a1.extract_span(spans[0]) == "MDNPERYMDMSGYQMDMQGRWMDAQGRFN")
        assert (a1.extract_span(spans[0]) == matching[0])
        assert (
            a1.extract_span(
                inverted[0]) == "NPFGQMWHGRQGHYPGYMSSHSMYGRNMYNPYHSHYASRH")
        assert (len(spans) == 5)
        assert (len(spans) == len(inverted) + 1)

    def test_equality_1(self):
        seq1 = sequence("ALKALI")
        seq2 = sequence("ALKALI")

        assert (seq1 == seq2)

    def test_equality_2(self):
        seq1 = sequence("ALKALI")
        seq2 = sequence("ALKALI")
        seq2.type = "rna"
        assert (seq1 != seq2)

    def test_lt_gt(self):
        seq1 = sequence("ALKALI")
        seq2 = sequence("BLKALI")
        seq3 = sequence("AALKALI")

        assert ((seq1 < seq2) == (seq1.seq < seq2.seq))
        assert ((seq3 < seq1) == (seq3.seq < seq1.seq))

    def test_calculate_weight(self):
        # TODO: test all amino acids?
        seq = sequence("ALKALI")
        assert (seq.calculate_weight() == 609.8099)

    def test_calculate_composition(self):
        seq = sequence("ALKALI")
        freqs = seq.calculate_composition()

        assert (freqs['A'] == 2)
        assert (freqs['K'] == 1)
        assert (freqs['W'] == 0)

    def test_find_open_reading_frames_1(self):
        seq = sequence("CGCTACGTCTTACGCTGGAGCTCTCATGGATCGGTTCGGTAGGGCTCGATCACATCGCTAGCCAT".replace("T",'U'))
        orfs = seq.find_open_reading_frames(min_orf_len=4)
        assert len(orfs) == 3

    def test_find_codon(self):
        seq = sequence("CGTUAGCGT")
        ind = seq.find_codon("UAG")
        assert ind == 3

        ind = seq.find_codon("AUG")
        assert ind == -1

        seq = sequence("ALKALI")
        ind = seq.find_codon("AUG")
        assert ind is None

    def test_find_index_after_alignment(self):
        seq = sequence("A-C--T-AGT")
        unaligned = sequence("TTTTTACTAGTGGGGG")
        
        idx = seq.find_index_after_alignment(unaligned, 6)

        assert idx == 2

    def test_slice(self):
        seq = sequence("AAAAGGGGTTTTAAAA")
        
        sub_seq = seq[4:8]
        assert sub_seq == 'GGGG'