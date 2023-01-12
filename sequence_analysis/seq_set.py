"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files.
"""
from Bio import SeqIO
from sequence_analysis.sequence import sequence
import pdb
from sequence_analysis.utils import aa_alphabet
from sequence_analysis.utils import dna_alphabet
from sequence_analysis.utils import rna_alphabet
from sequence_analysis.utils import diff_letters

class seq_set:
    def __init__(self, list_of_sequences=None, file_name=None, type=None):
        if list_of_sequences is None:
            list_of_sequences = []
        self.records = list_of_sequences
        self.type = type

        if len(list_of_sequences) == 0 and file_name is None:
            print("WARNING: sequence set initialized without sequences or file name to read from. Please use read_fasta() to read in sequences or use add_sequence() to add sequences.")

        if file_name is not None:
            self.read_fasta(file_name)

        if type is None and len(self.records) > 0:
            self.set_type()

    def write_fasta(self,file_name):
        # function to write all sequences within the list to a fasta file
        f = open(file_name,'w')
        for seq in self.records:
            s = seq.seq
            name = "> " + seq.name
            f.write(name)
            f.write('\n')
            for ct, char in enumerate(s):
                if ct % 79 == 0 and ct != 0:
                    f.write('\n')
                f.write(char)
            f.write('\n')
        f.close()

    def get_letters(self):
        all_letters = set()
        for seq in self.records:
            s = seq.seq.upper()
            all_letters = all_letters.union(set([*s]))
        return all_letters

    def set_type(self,type=None):
        # TODO: do further testing of type detection
        all_letters = self.get_letters()
        
        if type is not None:
            self.type = type
            # check that the assigned type is consistent
            return

        if len(diff_letters.intersection(all_letters)):
            self.type = 'protein'
        elif set(dna_alphabet).issubset(all_letters) or set(rna_alphabet).issubset(all_letters):
            if 'U' in all_letters:
                self.type = 'rna'
            else:
                # this could fail if we have RNA sequences that happens to have no Us, but that's unlikely
                self.type = 'dna'
        else:
            print("I can't assign a type, please specify manually.")


    def read_fasta(self,file_name):
        if len(self.records) > 0:
            print("Warning: overwriting existing data")
        for record in SeqIO.parse(file_name,"fasta"):
            seq = sequence(str(record.seq), record.name)
            self.records.append(seq)

        self.set_type()

        for seq in self.records:
            seq.type = self.type

    def add_sequence(self,sequences):
        if isinstance(sequences, list):
            for seq in sequences:
                assert(isinstance(seq, sequence))
                self.records.append(seq)
            old_type = self.type
            self.set_type()
            assert(old_type == self.type)
            self.records[-1].type = self.type
        elif isinstance(sequences, sequence):
            self.records.append(sequences)
            old_type = self.type
            assert(old_type == self.type)
            self.set_type()
            self.records[-1].type = self.type
        else:
            print("ERROR: incorrect format for adding sequence. Please provide a list of sequences or a sequence object. Sequence was not added.")


