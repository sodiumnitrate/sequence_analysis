from sequence_analysis import seq_set
import pdb

class TestClust:
    def test_simple_clustering(self):
        prots = seq_set.seq_set(file_name="aux_files/cluster_test.fasta")
        labels = prots.cluster(2,gap=0,gap_open=-15)
        assert(labels[0] == labels[1])
        assert(labels[1] == labels[2])
        assert(labels[3] == labels[4])
        assert(labels[4] == labels[5])

    def test_simple_clustering_2(self):
        prots = seq_set.seq_set(file_name="aux_files/cluster_test.fasta")
        prots.records[0], prots.records[-1] = prots.records[-1], prots.records[0]
        labels = prots.cluster(2,gap=0,gap_open=-15)

        assert(labels[0] == labels[3])
        assert(labels[3] == labels[4])
        assert(labels[1] == labels[2])
        assert(labels[2] == labels[5])