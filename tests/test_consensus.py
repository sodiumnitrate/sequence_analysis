from sequence_analysis.consensus import Consensus
from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence
import pdb

class TestConsensus:
    def test_init(self):
        seqs = seq_set(file_name="aux_files/cluster_test_2.fasta")
        con = Consensus(seqs)

        assert con.aligned_len is None
        assert not con.aligned 
        assert con.type == 'protein'
        assert con.n == 20

    def test_align(self):
        seqs = seq_set(file_name="aux_files/cluster_test_2.fasta")
        con = Consensus(seqs)
        con.align_seqs()

        assert con.aligned
        assert con.aligned_len is not None

    def test_calc_freq(self):
        seqs = seq_set(file_name="aux_files/cluster_test_2.fasta")
        con = Consensus(seqs)
        con.align_seqs()

        con.calculate_frequency()

        assert con.frequency is not None
        assert con.entropy is not None
        assert con.information is not None

        assert con.frequency.shape == (len(con.alphabet), len(con.sequences[0]))
        assert len(con.entropy) == len(con.sequences[0])

    def test_consensus(self):
        seqs = seq_set(file_name="aux_files/cluster_test_2.fasta")
        con = Consensus(seqs)
        con.align_seqs()

        con.calculate_frequency()

        con.get_consensus()
        assert con.consensus is not None

        con.get_consensus(non_gap_cutoff=33)
        assert con.consensus is not None
        assert len(con.consensus) == 61

    def test_get_codon_aa_consensus(self):
        seq1 = sequence("CAA")
        seq2 = sequence("GAA")
        seq3 = sequence("CAA")
        seq4 = sequence("CAT")
        seq5 = sequence("CAA")
        sset = seq_set([seq1, seq2, seq3, seq4, seq5])
        sset.type = "dna"

        con = Consensus(sset)
        con.get_consensus()
        possibilities = con.get_codon_aa_consensus(0)

        assert possibilities[0][0] == "Q"