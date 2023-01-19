from sequence_analysis.pairwise_alignment import pairwise_alignment

seq1 = "HEAGAWGHEE"
seq2 = "PAWHEAE"
new_alignment = pairwise_alignment(seq1,seq2,algorithm="smith-waterman")
new_alignment.align(verbose=True)

print(new_alignment.sequence1_aligned.seq)
print(new_alignment.sequence2_aligned.seq)
print(new_alignment.score)