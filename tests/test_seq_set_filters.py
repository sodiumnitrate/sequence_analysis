from sequence_analysis.seq_set import seq_set
import pdb

class TestFilters:
    def test_remove_duplicates(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        assert(len(prots.records)==6)
        prots.remove_duplicates()
        assert(len(prots.records)==3)

    def test_add_test(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        dnas = seq_set(file_name="aux_files/dna_ex.fasta")
        len_prots = prots.get_len()
        prots.add_set(dnas)
        assert(len_prots == prots.get_len())

    def test_filter_by_regex(self):
        prots = seq_set(file_name="aux_files/filter_test.fasta")
        prots.filter_by_pattern("MDM")
        assert(prots.get_len() == 2)

    def test_get_frequencies(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        unique_records, freqs = prots.get_frequencies()

        assert(len(unique_records) == 3)
        assert(freqs[0] == 2)
        assert(freqs[1] == 2)
        assert(freqs[2] == 2)

        prots.remove_duplicates()
        assert(len(unique_records) == prots.get_len())

    def test_get_frequencies_2(self):
        prots = seq_set(file_name="aux_files/dup_test_2.fasta")
        unique_records, freqs = prots.get_frequencies()

        assert(len(unique_records) == 3)
        assert(freqs[0] == 2)
        assert(freqs[1] == 2)
        assert(freqs[2] == 1)

    def test_sort(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        prots.alphabetize()
        assert(prots.records[0] == prots.records[1])

    def test_similarity_matrix(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        sim_matrix = prots.get_similarity_matrix()
        n = prots.get_len()
        m = (n-1) * n / 2
        assert(len(sim_matrix) == m)