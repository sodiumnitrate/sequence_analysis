"""
Basic tests for the seq_seq class.
"""
from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence

class TestSeqSet:
    def test_set_type(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')
        assert sequences.type == 'dna'

        sequences = seq_set(file_name="aux_files/pub.fasta")
        assert sequences.type == 'protein'

    def test_slice_1(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')

        seq = sequences[0]

        assert isinstance(seq, sequence)

        assert seq == sequences.records[0]

    def test_slice_2(self):
        sequences = seq_set(file_name='aux_files/dna_ex.fasta')

        new_set = sequences[:2]

        assert isinstance(new_set, seq_set)
        assert len(new_set) == 2

    def test_alphabetize(self):
        sequences = seq_set(file_name='aux_files/cluster_test.fasta')

        sequences.alphabetize()
        assert sequences[0].name == 'seq6'

    def test_unique_naming(self):
        sequences = seq_set(file_name='aux_files/dup_test_2.fasta')

        sequences.remove_duplicates()
        assert len(sequences) == 3
        assert sequences[0].name == "seq2_seq5"

    def test_add_gaps_to_align_pattern(self):
        seq1 = sequence("MMMMMYYYMMMMM")
        seq2 = sequence("MMYYYMMMMMMMMMMM")
        n = len(seq2)
        seqs = seq_set(list_of_sequences=[seq1, seq2])
        seqs.add_gaps_to_align_pattern("YYY")

        assert len(seqs[1]) == n + 3
        assert seqs[1][0] == '-'

    def test_seq_logo_dna(self):
        seq1 = sequence("AGCTAGCT")
        seq2 = sequence("AGCTACGT")
        seq3 = sequence("TAGTCGCT")
        seq4 = sequence("TAGCCGCA")

        seqs = seq_set(list_of_sequences=[seq1, seq2, seq3, seq4])

        alpha, height = seqs.create_sequence_logo((0,8))

        assert len(alpha) == 4
        assert len(height[0]) == 8

        assert height[0][0] + height[1][0] + height[2][0] + height[3][0] == height[0][1] + height[1][1] + height[2][1] + height[3][1]


    def test_seq_logo_rna(self):
        seq1 = sequence("AGCUAGCU")
        seq2 = sequence("AGCUACGU")
        seq3 = sequence("UAGUCGCU")
        seq4 = sequence("UAGCCGCA")

        seqs = seq_set(list_of_sequences=[seq1, seq2, seq3, seq4])

        alpha, height = seqs.create_sequence_logo((0,8))

        assert len(alpha) == 4
        assert len(height[0]) == 8

        assert height[0][0] + height[1][0] + height[2][0] + height[3][0] == height[0][1] + height[1][1] + height[2][1] + height[3][1]

    def test_seq_logo_prot(self):
        seq1 = sequence("MDMALKALI")
        seq2 = sequence("ALKALIMDM")
        seq3 = sequence("MALKALIDM")
        seq4 = sequence("AL-DMKALI")

        seqs = seq_set(list_of_sequences=[seq1, seq2, seq3, seq4])

        alpha, height = seqs.create_sequence_logo((0,8))

        assert len(alpha) == 20
        assert len(height[0]) == 8

        height_sums = [0 for _ in range(8)]
        for i in range(8):
            for j in range(20):
                height_sums[i] += height[j][i]

        # I expect the same information content in position 0 and position 7
        assert height_sums[0] == height_sums[-1]