from sequence_analysis import seq_set

dna = seq_set.seq_set(file_name="example_files/sample_dna.fasta")

seq = dna.records[0]

true_sequence = seq.six_frame_check("MDM")
