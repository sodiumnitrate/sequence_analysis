"""
Defines the Consensus class. Holds a list of sequences. Aligns them if need be
and finds the consensus. 

TODO:
- Can find A.A. consensus from DNA/RNA seqs.
"""

import numpy as np

from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence
from sequence_analysis.utils import aa_alphabet, rna_alphabet, dna_alphabet
from sequence_analysis.msa import MSA

class Consensus:
    """
    Consensus object to align (if needed) a set of sequences and determine
    consensus.
    """
    def __init__(self, sequences, aligned=False):
        self.sequences = sequences
        self.aligned = aligned
        self.aligned_len = None
        self.check_data()
        self.type = self.sequences.type
        self.n = None
        self.alphabet = None
        self.check_type()

        self.frequency = None
        self.information = None
        self.entropy = None
        self.total_frequency = None

        self.consensus = None

    def check_data(self):
        """
        Make sure that the input data makes sense.
        """
        if not isinstance(self.sequences, seq_set):
            print("ERROR: sequences must be of type seq_set.")
            raise TypeError
        if self.aligned:
            lens = set([len(a) for a in self.sequences])
            if len(lens) != 1:
                print("ERROR: aligned sequences can't have different lengths.")
                raise ValueError
            else:
                self.aligned_len = list(lens)[0]

    def align_seqs(self):
        """
        Align sequences using the MSA object.
        """
        msa = MSA(self.sequences)
        msa.align(clean=True)
        self.sequences = msa.aligned_sequences
        self.aligned_len = len(self.sequences[0])
        self.aligned = True

    def check_type(self):
        """
        Check type of sequences and assign alphabet.
        """
        if self.type == "protein":
            self.n = 20
            self.alphabet = aa_alphabet + ["X"]
        elif self.type == "rna":
            self.n = 4
            self.alphabet = rna_alphabet + ["N"]
        elif self.type == "dna":
            self.n = 4
            self.alphabet = dna_alphabet + ["N"]
        else:
            print("ERROR: sequence type not detected.")
            raise TypeError
        self.alphabet = {a:i for i,a in enumerate(self.alphabet)}

    def calculate_frequency(self):
        """
        Calculate position-specific frequences, entropy, and information content.
        """
        if not self.aligned:
            self.align_seqs()

        # we are using ones so that there are no zero frequencies in the end
        # see Lawrence et al. Science 1993
        freqs = np.ones((len(self.alphabet), self.aligned_len))
        for seq in self.sequences:
            for i, char in enumerate(seq.seq):
                if char != '-':
                    idx = self.alphabet[char]
                    freqs[idx,i] += 1

        tot_freq = np.sum(freqs, axis=0)
        rel_freq = freqs / tot_freq
        e_ni = 1./ np.log(2) * (self.n - 1) / (2 * tot_freq)
        H_i = -1 * np.sum(rel_freq * np.log2(rel_freq), axis=0)
        R_i = np.log2(self.n) - (H_i + e_ni)

        self.frequency = rel_freq
        self.entropy = H_i
        self.information = R_i
        self.total_frequency = tot_freq

    def get_consensus(self, non_gap_cutoff=1):
        """
        Get consensus, based on most likely at each position.
        """
        rev_alpha = {val:key for key,val in self.alphabet.items()}
        if self.frequency is None or self.entropy is None or self.information is None:
            self.calculate_frequency()

        consensus = ""
        for i in range(self.aligned_len):
            if self.total_frequency[i] >= non_gap_cutoff:
                idx = np.argmax(self.frequency[:,i])
                char = rev_alpha[idx]
                consensus += char

        self.consensus = sequence(consensus)

    def get_codon_aa_consensus(self, idx):
        """
        Given an index, get codon starting at that index. Then, get probabilities of all
        64 combinations, with corresponding amino acid.
        
        Skips char=N.
        """
        if idx+2 >= self.aligned_len:
            print(f"ERROR: a codon starting at index={idx} out of range.")
            raise IndexError

        if self.type == "protein":
            print("ERROR: protein sequences don't have codons.")
            raise TypeError

        if self.frequency is None:
            print("WARNING: consensus not found -- proceeding with default options.")
            self.get_consensus()

        rev_alpha = {val:key for key,val in self.alphabet.items()}

        pos1 = self.frequency[:,idx]
        pos2 = self.frequency[:,idx+1]
        pos3 = self.frequency[:,idx+2]

        possibilities = {}
        for i1 in range(4):
            for i2 in range(4):
                for i3 in range(4):
                    prob = pos1[i1]*pos2[i2]*pos3[i3]
                    codon = sequence(rev_alpha[i1] + rev_alpha[i2] + rev_alpha[i3])
                    codon.type = self.type
                    aa = codon.translate().seq
                    if aa in possibilities:
                        possibilities[aa] += prob
                    else:
                        possibilities[aa] = prob

        return sorted(possibilities.items(), key=lambda x:x[1], reverse=True)

    def get_protein_consensus_from_nucleotide(self, frame=0, strand=1):
        """
        Given a dna or rna sequence, get consensus to form a protein consensus.
        """
        if self.frequency is None or self.entropy is None or self.information is None:
            self.calculate_frequency()

        idx = frame
        while idx + 2 < self.aligned_len:
            pos = self.get_codon_aa_consensus(idx)
            