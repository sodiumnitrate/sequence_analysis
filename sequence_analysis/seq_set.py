"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files, and various filtering methods.
"""
import re
import copy
import numpy as np
from Bio import SeqIO
from sequence_analysis.utils import add_dicts
from sequence_analysis.sequence import sequence
from sequence_analysis.utils import rna_alphabet, diff_letters
from sequence_analysis.utils import write_in_columns, aa_alphabet, dna_alphabet
from sequence_analysis.utils import merge_dicts, fasta_or_phylip
import fasta_reader_cpp

class seq_set:
    def write_phylip(self, file_name, mode='sequential'):
        """
        Write sequences into a phylip file.

        Currently only supporting sequential mode (which codeml can process).
        """
        file = open(file_name, 'w', encoding="utf-8")
        n_seqs = len(self)
        n_chars = list(set([len(seq) for seq in self.records]))
        if len(n_chars) != 1:
            print("ERROR: sequences must be of the same length.")
            raise ValueError

        n_chars = n_chars[0]

        if mode == "sequential":
            file.write(f"{n_seqs}\t{n_chars}\n")
        elif mode == "interleaved":
            file.write(f"{n_seqs}\t{n_chars}\tI\n")

        if mode == "sequential":
            for seq in self.records:
                file.write(f"{seq.name}    {seq.seq}\n")
        elif mode == "interleaved":
            # write in 6 blocks of 10
            n_name_chars = max([len(s.name) for s in self])
            total_lines = int(np.ceil(n_chars / 60))
            for i in range(total_lines):
                for seq in self.records:
                    if i == 0:
                        file.write(f"{seq.name}{' ' * (n_name_chars - len(seq.name) + 2)}")
                    else:
                        file.write(" " * (n_name_chars + 2))
                    seq_string = " ".join([seq[i*60+j*10:min(i*60+(j+1)*10, len(seq))] for j in range(6)])
                    file.write(seq_string)
                    file.write("\n")
                file.write("\n")
        file.close()



    def read_fastq(self, file_name):
        """
        Function to read sequences into seq_set from a .fastq file.
        """
        seq_str = ""
        quality_str = ""
        prev_char = None
        with open(file_name, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    if prev_char is not None:
                        seq = sequence(seq_str, seq_name)
                        seq.quality = quality_str
                        self.records.append(seq)
                        seq_str = ""
                        quality_str = ""
                    seq_name = line.strip()[1:]
                    prev_char = '@'
                elif line.startswith('+'):
                    prev_char = '+'
                else:
                    if prev_char == '@':
                        seq_str = seq_str + line.strip()
                    else:
                        quality_str = quality_str + line.strip()
        seq = sequence(seq_str, seq_name)
        seq.quality = quality_str
        self.records.append(seq)
        return

    def write_fastq(self, file_name):
        """
        Write sequences into fastq file.
        """
        f = open(file_name, 'w')
        for seq in self.records:
            f.write(f"@{seq.name}\n")
            f.write(f"{seq.seq}\n")
            f.write("+\n")
            f.write(f"{seq.quality}\n")

        f.close()

    def read_phylip(self, file_name):
        """
        Function to read sequences from phylip files.
        """
        self.records = []
        with open(file_name, 'r') as f:
            header = f.readline()
            n_seqs = int(header.split()[0])
            n_chars = int(header.split()[1])
            if 'I' in header:
                # interleaved
                seq_ct = 0
                names = []
                seqs = ["" for i in range(n_seqs)]
                for line in f:
                    if line.strip():
                        ls = line.split()
                        if seq_ct < n_seqs:
                            names.append(ls[0])
                            seq_portion = "".join(ls[1:])
                            seqs[seq_ct] = seq_portion
                        else:
                            i = seq_ct % n_seqs
                            seqs[i] = seqs[i] + "".join(ls)
                    seq_ct += 1

                for i in range(n_seqs):
                    seq = sequence(seqs[i], names[i])
                    self.records.append(seq)

            else:
                # sequential
                for line in f:
                    if line.strip():
                        ls = line.split()
                        if len(ls[1]) != n_chars:
                            print("ERROR: sequence length does not match the first line of the file.")
                            raise ValueError
                        seq = sequence(ls[1], ls[0])
                        self.records.append(seq)

        if len(self.records) != n_seqs:
            print("ERROR: number of sequences does not match the first line of the file.")
            raise ValueError

        self.set_type()
        for seq in self.records:
            seq.type = self.type