"""
Unit tests for the SamEntry class.
"""
from sequence_analysis.sam_entry import SamEntry

import pdb

class TestSamEntry:
    def test_init(self):
        sam = SamEntry("SRR18071790.63502	0	CM042177.1	22464269	255	12M1D61M	*	0	0	CCCACTTGCCACTGTGAGCGTTCCTCAACCAATAAACGCGGTCAGAGCTCAATTCAGAGAAGATGAAATTGGT	<AAFFFFJFJJ<JFFJ-FJJFAJJJJ<AFJJJF--AF-F-7AJA<AF-<F7FFJ<J7A-FJ-F7AJJJ<-777	NH:i:1	HI:i:1	AS:i:65	nM:i:1")

        assert sam.q_name == "SRR18071790.63502"
        assert sam.bitwise_flag == 0
        assert sam.r_name == "CM042177.1"
        assert sam.pos == 22464268
        assert sam.mapping_quality == 255
        assert sam.cigar_string == "12M1D61M"
        assert sam.r_next == "*"
        assert sam.p_next == "0"
        assert sam.t_len == "0"
        assert sam.sequence.seq == sam.sequence_string
        assert len(sam.quality_33) == len(sam.sequence_string)
        assert len(sam.probability_of_incorrect_read) == len(sam.sequence_string)
        assert 'nM:i:1' in sam.other_flags

    def test_init_reverse(self):
        sam = SamEntry("SRR18071790.63263	272	CM042177.1	22898393	3	100M	*	0	0	CGGAAGGTCCATGTTCAACCAGGGACACAGCATGGACAGTCAACGCTACGGCGGCTGGATGGACAACCCCGAGAGGTACATGGACATGTCTGGCTACCAG	JJJJJJJ<<JFJJJJFJJJFAF77J-<FJJJJJA-F-FF7FFFJJJJJFAJ<JF7AF7F<JFJJJJFJJJF7JFFJJFFA7JJJJJJJAJJJJJAFFFAA	NH:i:2	HI:i:2	AS:i:82	nM:i:8")

        assert sam.sequence.reverse_complement().seq == sam.sequence_string

