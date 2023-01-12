from sequence_analysis import seq_set
from sequence_analysis import sequence

class TestIO:
    def test_read_fasta(self):
        prot_test = seq_set.seq_set(file_name="aux_files/test.fasta")
        assert(isinstance(prot_test, seq_set.seq_set))
        assert(isinstance(prot_test.records, list))
        assert(len(prot_test.records) == 3)
        assert(prot_test.type == 'protein')

    def test_add_sequence(self):
        prot_test = seq_set.seq_set()
        new_sequence = sequence.sequence("ALKALI",name="test_add")
        prot_test.add_sequence(new_sequence)
        assert(new_sequence.type == 'protein')
        assert(prot_test.type == 'protein')
        assert(len(prot_test.records) == 1)

    def test_read_dna_fasta(self):
        dna_test = seq_set.seq_set(file_name="aux_files/dna_ex.fasta")
        assert(dna_test.type == 'dna')
        assert(len(dna_test.records) == 2)