from sequence_analysis import seq_set

class TestClust:
    def test_simple_clustering(self):
        prots = seq_set.seq_set(file_name="aux_files/cluster_test.fasta")
        labels = prots.cluster(2,gap=0,gap_open=-15)
        assert(labels[0] == 0)
        assert(labels[1] == 0)
        assert(labels[2] == 0)
        assert(labels[3] == 1)
        assert(labels[4] == 1)
        assert(labels[5] == 1)