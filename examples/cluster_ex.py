from sequence_analysis.seq_set import seq_set
import matplotlib.pyplot as plt
import pdb

prots = seq_set(file_name="../tests/aux_files/cluster_test_2.fasta")
prots.cluster(3)
prots.visualize_clusters()
plt.show()
