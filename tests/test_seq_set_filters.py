from sequence_analysis.seq_set import seq_set

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
        assert(unique_records.get_len() == 3)
        assert(freqs[0] == 2)
        assert(freqs[1] == 2)
        assert(freqs[2] == 2)
