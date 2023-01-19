from sequence_analysis.seq_set import seq_set

class TestFilters:
    def test_remove_duplicates(self):
        prots = seq_set(file_name="aux_files/dup_test.fasta")
        assert(len(prots.records)==6)
        prots.remove_duplicates()
        assert(len(prots.records)==3)