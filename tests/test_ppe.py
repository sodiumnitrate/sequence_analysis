"""
Unit tests for the PPE class and its methods.
"""
from sequence_analysis.seq_set import seq_set
from sequence_analysis.ppe import PPE

import pdb

class TestPPE:
    def test_ppe_init_seq_set(self):
        sequences = seq_set(file_name='aux_files/test_2.fasta')
        motif_finder = PPE(sequences)

        assert isinstance(motif_finder.sequences, list)
        assert all([isinstance(s, list) for s in motif_finder.sequences])
        assert all([all([isinstance(s, str) for s in seq]) for seq in motif_finder.sequences])

    def test_ppe_init_list_of_seqs(self):
        sequences = seq_set(file_name='aux_files/test_2.fasta')
        motif_finder = PPE(sequences.records)

        assert isinstance(motif_finder.sequences, list)
        assert all([isinstance(s, list) for s in motif_finder.sequences])
        assert all([all([isinstance(s, str) for s in seq]) for seq in motif_finder.sequences])

    def test_ppe_init_list_of_strings(self):
        sequences = seq_set(file_name='aux_files/test_2.fasta')
        sequences = [seq.seq for seq in sequences]

        motif_finder = PPE(sequences)

        assert isinstance(motif_finder.sequences, list)
        assert all([isinstance(s, list) for s in motif_finder.sequences])
        assert all([all([isinstance(s, str) for s in seq]) for seq in motif_finder.sequences])

    def test_ppe_init_list_of_list_of_strings(self):
        sequences = seq_set(file_name='aux_files/test_2.fasta')
        sequences = [list(seq.seq) for seq in sequences]

        motif_finder = PPE(sequences)

        assert isinstance(motif_finder.sequences, list)
        assert all([isinstance(s, list) for s in motif_finder.sequences])
        assert all([all([isinstance(s, str) for s in seq]) for seq in motif_finder.sequences])

    def test_generate_alphabet_1(self):
        sequences = ['AAACAAACALKAA','LKAAAACLKAAAAC']

        motif_finder = PPE(sequences)

        assert len(motif_finder.alphabet) == 4

    def test_generate_alphabet_2(self):
        sequences = [['AA','ACAAACALKAA'],['LKAA','AACLKAAAAC']]

        motif_finder = PPE(sequences)

        assert len(motif_finder.alphabet) == 4
        assert all([all([j in motif_finder.alphabet for j in s]) for s in sequences])

    def test_find_most_common_twomer(self):
        sequences = ['AAACAAACALKAA','LKAAAACLKAAAAC']
        motif_finder = PPE(sequences)

        motif_finder.find_most_common_twomer()
        assert motif_finder.merge_operations[-1] == ('A', 'A')

    def test_merge_once(self):
        sequences = ['AAACAAACALKAA','LKAAAACLKAAAAC']
        motif_finder = PPE(sequences)

        motif_finder.find_most_common_twomer()
        motif_finder.merge_and_update_alphabet()

        assert 'AA' in motif_finder.alphabet

    def test_iterate(self):
        sequences = ['AAACAAACALKAA','LKAAAACLKAAAAC']
        motif_finder = PPE(sequences)

        motif_finder.iterate()

        assert motif_finder.n_iter > 1
        assert motif_finder.frequency_of_seqs_with_most_common >= motif_finder.frequency_of_most_common_motif

        assert motif_finder.motifs[0][0] == 'LKAAAAC'
        assert motif_finder.motifs[0][1] == 2

    def test_gen_vec(self):
        sequences = ['AAACAAACALKAA','LKAAAACLKAAAAC']
        motif_finder = PPE(sequences)

        motif_finder.iterate()

        motif_finder.generate_vectors()

        assert len(motif_finder.vectors) == len(sequences)
        assert len(motif_finder.vectors[0]) == len(motif_finder.alphabet)