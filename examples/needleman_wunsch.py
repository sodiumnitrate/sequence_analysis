from sequence_analysis.pairwise_alignment import pairwise_alignment

seq1 = "ALKALIMDMKALI"
seq2 = "ALKAKIWMKAL"
new_alignment = pairwise_alignment(seq1,seq2)
new_alignment.align()

print(new_alignment.sequence1_aligned.seq)
print(new_alignment.sequence2_aligned.seq)
print(new_alignment.score)