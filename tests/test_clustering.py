"""
Script that tests various methods of the clustering class.
"""

from sequence_analysis import seq_set


class TestClust:
    """
    Class that contains various tests.
    """

    def test_simple_clustering(self):
        """Function that tests simple clustering with 6 sequences."""
        prots = seq_set.seq_set(file_name="aux_files/cluster_test.fasta")
        prots.cluster(2, gap=0, gap_open=-15)

        labels = prots.cluster_labels

        assert (labels[0] == labels[1])
        assert (labels[1] == labels[2])
        assert (labels[3] == labels[4])
        assert (labels[4] == labels[5])

    def test_simple_clustering_2(self):
        """Function that tests simple clustering with 6."""
        prots = seq_set.seq_set(file_name="aux_files/cluster_test.fasta")
        prots.records[0], prots.records[-1] = prots.records[-1], prots.records[0]
        prots.cluster(2, gap=0, gap_open=-15)

        labels = prots.cluster_labels

        assert (labels[0] == labels[3])
        assert (labels[3] == labels[4])
        assert (labels[1] == labels[2])
        assert (labels[2] == labels[5])
