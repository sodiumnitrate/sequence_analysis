"""
Tests for the seq_seq class.
"""

from sequence_analysis.seq_set import seq_set

class TestFilters:
    def test_len(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        assert (len(prots) == 6)

    def test_remove_duplicates(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        assert (len(prots.records) == 6)
        prots.remove_duplicates()
        assert (len(prots.records) == 3)

    def test_add_test(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        dnas = seq_set(file_name="aux_files/dna_ex.fasta")
        len_prots = prots.get_len()
        prots.add_set(dnas)
        assert (len_prots == prots.get_len())

    def test_filter_by_regex(self):
        prots = seq_set(file_name="aux_files/filter_test.fasta")
        prots.filter_by_pattern("MDM")
        assert (prots.get_len() == 2)

    def test_get_frequencies(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        unique_records, freqs = prots.get_frequencies()

        assert (len(unique_records) == 3)
        assert (freqs[0] == 2)
        assert (freqs[1] == 2)
        assert (freqs[2] == 2)

        prots.remove_duplicates()
        assert (len(unique_records) == prots.get_len())

    def test_get_frequencies_2(self):
        prots = seq_set(file_name="aux_files/dup_test_2.fasta")
        unique_records, freqs = prots.get_frequencies()

        assert (len(unique_records) == 3)
        assert (freqs[0] == 2)
        assert (freqs[1] == 2)
        assert (freqs[2] == 1)

    def test_sort(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        prots.alphabetize()
        assert (prots.records[0] == prots.records[1])

    def test_similarity_matrix(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        prots.get_similarity_matrix()
        sim_matrix = prots.sim_matrix
        n = prots.get_len()
        m = (n - 1) * n / 2
        assert (len(sim_matrix) == m)

    def test_filter_by_frequency(self):
        prots = seq_set(file_name="aux_files/dup_test_2.fasta")
        prots.filter_by_frequency(threshold=2)
        assert (prots.get_len() == 2)

    def test_calculate_composition_dna_collate(self):
        dnas = seq_set(file_name="aux_files/dna_ex.fasta")
        freqs = dnas.calculate_composition(collate=True)

        assert (freqs['A'] == 4)
        assert (freqs['C'] == 1)
        assert (freqs['T'] == 4)
        assert (freqs['G'] == 2)

        assert ('W' not in freqs.keys())

    def test_calculate_composition_dna_individual(self):
        dnas = seq_set(file_name="aux_files/dna_ex.fasta")
        freqs = dnas.calculate_composition()

        assert (isinstance(freqs, type([])))
        assert (len(freqs) == dnas.get_len())

        assert (freqs[0]['A'] == 2)

    def test_calculate_composition_protein_collate(self):
        prots = seq_set(file_name="aux_files/filter_test.fasta")
        freqs = prots.calculate_composition(collate=True)

        assert (freqs['W'] == 19)

    def test_filter_by_weight(self):
        prots = seq_set(file_name="aux_files/test.fasta")
        prots.filter_by_weight(threshold=7000)

        assert (prots.get_len() == 2)

    def test_filter_by_six_frame(self):
        rnas = seq_set(file_name="aux_files/frame_shift.fasta")
        prots = rnas.filter_by_six_frame_check_pattern("MDM")
        assert (len(prots) == 3)