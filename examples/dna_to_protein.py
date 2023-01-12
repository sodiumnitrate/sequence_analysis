from sequence_analysis import seq_set
from sequence_analysis import sequence
import pdb

dna = seq_set.seq_set(file_name="example_files/sample_dna.fasta")

seq = dna.records[0]

selected, true_sequence = seq.six_frame_check("MDM")