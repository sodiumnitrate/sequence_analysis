from sequence_analysis import seq_set
from sequence_analysis import sequence

# create sequence set object
test = seq_set.seq_set()

test.read_fasta("example_files/test.fasta")

# number of sequences
print(f"Number of sequences: {len(test.records)}")
print(f"Type of sequences: {test.type}")
# create new sequence
new_sequence = sequence.sequence("ALKALI", name="test_add")
test.add_sequence(new_sequence)

print(f"Sequence added. Now we have: {len(test.records)}")
print(f"The added sequence is '{test.records[-1].seq}'")

# write the new set to fasta file
test.write_fasta("example_files/new_test.fasta")

# load a .fasta file containing DNA sequences and write out the type
dna_test = seq_set.seq_set(file_name="example_files/dna_ex.fasta")
print(
    f"Loaded DNA fasta file. There are {len(dna_test.records)} sequences of type {dna_test.type}")
