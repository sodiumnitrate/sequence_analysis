from sequence_analysis import seq_set
from sequence_analysis import sequence


class TestIO:
    def test_read_fasta(self):
        prot_test = seq_set.seq_set(file_name="aux_files/test.fasta")
        assert (isinstance(prot_test, seq_set.seq_set))
        assert (isinstance(prot_test.records, list))
        assert (len(prot_test.records) == 3)
        assert (prot_test.type == 'protein')

    def test_add_sequence(self):
        prot_test = seq_set.seq_set()
        new_sequence = sequence.sequence("ALKALI", name="test_add")
        prot_test.add_sequence(new_sequence)
        assert (new_sequence.type == 'protein')
        assert (prot_test.type == 'protein')
        assert (len(prot_test.records) == 1)

    def test_add_sequence_by_list(self):
        prot_test = seq_set.seq_set()
        new_sequence = sequence.sequence("ALKALI", name="test_add")
        new_sequence2 = sequence.sequence("ALKALIMDM", name="test_add2")
        prot_test.add_sequence([new_sequence, new_sequence2])
        assert (new_sequence.type == "protein")
        assert (new_sequence2.type == "protein")
        assert (prot_test.get_len() == 2)

    def test_read_dna_fasta(self):
        dna_test = seq_set.seq_set(file_name="aux_files/dna_ex.fasta")
        assert (dna_test.type == 'dna')
        assert (len(dna_test.records) == 2)

    def test_write_fasta(self):
        prot_test = seq_set.seq_set(file_name="aux_files/test.fasta")
        prot_test.write_fasta("aux_files/write_test.fasta")
        test_read = seq_set.seq_set(file_name="aux_files/write_test.fasta")

        # make sure the two files are the same
        assert (prot_test.get_len() == test_read.get_len())

        for i in range(prot_test.get_len()):
            s1 = prot_test.records[i].seq
            s2 = test_read.records[i].seq
            assert (s1 == s2)
