"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files, and various filtering methods.
"""
from Bio import SeqIO
from sequence_analysis.sequence import sequence
from sequence_analysis.utils import dna_alphabet
from sequence_analysis.utils import rna_alphabet
from sequence_analysis.utils import diff_letters
from sequence_analysis.pairwise_alignment import pairwise_alignment
import time
from sklearn.cluster import SpectralClustering
import numpy as np
from sequence_analysis.utils import add_dicts

class seq_set:
    """This class holds a list of sequence objects of a given type."""
    def __init__(self, list_of_sequences=None, file_name=None, type=None):
        if list_of_sequences is None:
            list_of_sequences = []
        self.records = list_of_sequences
        self.type = type

        self.sim_matrix = None

        if len(list_of_sequences) == 0 and file_name is None:
            print("WARNING: sequence set initialized without sequences or file name to read from. Please use read_fasta() to read in sequences or use add_sequence() to add sequences.")

        if file_name is not None:
            self.read_fasta(file_name)

        if type is None and len(self.records) > 0:
            self.set_type()

    def __str__(self):
        """__str__ function for sequence set (seq_set) object."""
        return f"Sequence set object with {len(self)} sequences of type {self.type}"

    def __repr__(self):
        """"__repr__ function for sequence set (seq_set) object."""
        return f"<seq_set object of size {len(self)} and type {self.type} at {hex(id(self))}>"

    def __len__(self):
        """Overwrites __len__ to return number of sequences within set."""
        return len(self.records)
    
    def get_len(self):
        """Function that returns the number of sequences within seq_set."""
        # TODO: switch to len() in all instances and get rid of this function.
        return len(self.records)

    def write_fasta(self,file_name):
        """"Function to write all sequences within seq_set to a fasta file."""
        f = open(file_name,'w')
        for seq in self.records:
            s = seq.seq
            name = ">" + seq.name
            f.write(name)
            f.write('\n')
            for ct, char in enumerate(s):
                if ct % 79 == 0 and ct != 0:
                    f.write('\n')
                f.write(char)
            f.write('\n')
        f.close()

    def filter_by_six_frame_check_pattern(self, regex, overwrite_frame_shifted=True):
        """
        Function that checks pattern in the translated sequence and returns
        the translated protein sequence.

        If overwrite_frame_shifted is True, sequences in the set will be 
        overwritten by the correct frame-shifted rna or dna sequence.
        """
        assert(self.type == 'rna' or self.type == 'dna')

        true_records = []
        untranslated_records = []
        for seq in self.records:
            true, untranslated = seq.six_frame_check(regex)
            if true is not None:
                true_records.append(true)
                untranslated_records.append(untranslated)

        if overwrite_frame_shifted:
            self.records = untranslated_records

        return true_records

    def remove_duplicates(self,alphabetize=True):
        """
        Function that removes duplicates within seq_set.
        Also alphabetizes records. Keeps the first occurence in the set.
        
        NOTE: this is allowed for by overwriting the eq operator in the sequence class.
        """
        if not alphabetize:
            new_records = []
            for seq in self.records:
                if seq not in new_records:
                    new_records.append(seq)

            self.records = new_records
        else:
            self.alphabetize()
            unique_records = []
            for i in range(len(self.records)):
                if i == 0:
                    unique_records.append(self.records[0])
                else:
                    if self.records[i] != self.records[i-1]:
                        unique_records.append(self.records[i])

            self.records = unique_records

    def get_frequencies(self):
        """Function that gets frequencies of all unique sequences."""
        records = self.records
        records.sort()

        unique_records = []
        frequencies = []
        for i in range(len(records)):
            if i == 0:
                count = 1
                unique_records.append(records[0])
            else:
                if records[i] == records[i-1]:
                    count += 1
                else:
                    frequencies.append(count)
                    count = 1
                    unique_records.append(records[i])

        frequencies.append(count)
        
        return unique_records, frequencies

    def alphabetize(self):
        """Function to sort sequences within seq_set in alphabetical order."""
        records = self.records
        records.sort()
        self.records = records

    def get_letters(self):
        """"Function to get all unique letters in the seq_set sequences."""
        all_letters = set()
        for seq in self.records:
            s = seq.seq.upper()
            all_letters = all_letters.union(set([*s]))
        return all_letters

    def set_type(self,type=None):
        """Function to detect and set the type of sequences in seq_set."""
        # TODO: do further testing of type detection
        all_letters = self.get_letters()
        
        if type is not None:
            self.type = type
            # TODO: check that the assigned type is consistent
            return

        if len(diff_letters.intersection(all_letters)):
            self.type = 'protein'
        elif set(dna_alphabet).issubset(all_letters) or set(rna_alphabet).issubset(all_letters):
            if 'U' in all_letters:
                self.type = 'rna'
            else:
                # TODO: this could fail if we have RNA sequences that happens to have no Us, but that's unlikely
                self.type = 'dna'
        else:
            print("I can't assign a type, please specify manually.")


    def read_fasta(self,file_name):
        """Function to read sequences into seq_set from a .fasta file."""
        if len(self.records) > 0:
            print("Warning: overwriting existing data")
        for record in SeqIO.parse(file_name,"fasta"):
            seq = sequence(str(record.seq), record.name)
            self.records.append(seq)

        self.set_type()

        for seq in self.records:
            seq.type = self.type

    def add_sequence(self,sequences):
        """Function that adds another sequence to seq_set."""
        # TODO: refactor, because you have an add_set method as well
        if isinstance(sequences, list):
            for seq in sequences:
                assert(isinstance(seq, sequence))
                self.records.append(seq)
            old_type = self.type
            self.set_type()
            if old_type != None:
                assert(old_type == self.type)
            self.records[-1].type = self.type
        elif isinstance(sequences, sequence):
            self.records.append(sequences)
            old_type = self.type
            self.set_type()
            if old_type != None:
                assert(old_type == self.type)
            self.records[-1].type = self.type
        else:
            print("ERROR: incorrect format for adding sequence. Please provide a list of sequences or a sequence object. Sequence was not added.")

    def add_set(self,set2):
        """Function that merges two sets of sequences (two seq_sets)."""
        # make sure the set being added is of type seq_set
        if not isinstance(set2,seq_set):
            print("ERROR: cannot add the two sets. Set2 is not of type seq_set.")
            return

        # make sure sequence types are the same
        if set2.type != self.type:
            print("ERROR: cannot add two sets. They are not of the same type.")
            return

        # add set2's records to the list of records 
        self.records += set2.records

    def filter_by_pattern(self,regex):
        """Filter sequences in seq_set by the given regex pattern."""
        new_records = []
        for s in self.records:
            if s.check_for_pattern(regex):
                new_records.append(s)

        self.records = new_records

    def get_similarity_matrix(self,algorithm="biopython-global",use_blosum_50=False,match=2,unmatch=-1,gap=-0.1,gap_open=-0.5,verbose=False):
        """
        Function that calculates similarity matrix for the sequences in seq_set.

        For a seq_set containing N sequences, calculates (N-1)N/2 pairwise alignment scores.
        """
    
        # pairwise alignment between all pairs of sequences
        n = self.get_len()
        similarity_matrix = {}
        for i in range(n):
            for j in range(i+1,n):
                start = time.time()
                alignment = pairwise_alignment(self.records[i],self.records[j],algorithm=algorithm,use_blosum_50=use_blosum_50,match=match,unmatch=unmatch,gap=gap,gap_open=gap_open)
                alignment.align()
                similarity_matrix[(i,j)] = alignment.score
                end = time.time()

                if verbose:
                    print(i,j,end-start)

        self.sim_matrix = similarity_matrix

    def filter_by_frequency(self,threshold=10):
        """
        Function that deletes sequences that have a frequency lower than the 
        given threshold.

        By construction, does not return duplicates of the sequences.
        """
        recs, frequencies = self.get_frequencies()
        new_recs = []
        for i, rec in enumerate(recs):
            if frequencies[i] >= threshold:
                new_recs.append(rec)

        self.records = new_recs

    def remove_before_pattern(self, regex, verbose=False):
        """
        Function that, for each sequence in seq_set, removes the letters before the 
        matching regex pattern.
        """
        for seq in self.records:
            seq.remove_before_pattern(regex,verbose=verbose)


    def cluster(self,n_clusters,algorithm="biopython-global", use_blosum_50=False, match=2, unmatch=-1, gap=-0.1, gap_open=-0.5, verbose=False):
        """
        Function that clusters the sequences in seq_set.

        Similarity matrix is calculated based on pairwise alignments, and then used to 
        perform spectral clustering.
        """
        if self.sim_matrix is None:
            self.get_similarity_matrix(algorithm=algorithm, use_blosum_50=use_blosum_50, match=match, unmatch=unmatch, gap=gap, gap_open=gap_open, verbose=verbose)
        else:
            n = self.get_len()
            N = len(self.sim_matrix)
            if N != (n-1)*n/2:
                self.get_similarity_matrix(algorithm=algorithm, use_blosum_50=use_blosum_50, match=match, unmatch=unmatch, gap=gap, gap_open=gap_open, verbose=verbose)

        sim = self.sim_matrix

        n = self.get_len()
        sim_matrix = np.zeros((n,n))
        for el in sim.keys():
            sim_matrix[el[0],el[1]] = sim[el]
            sim_matrix[el[1],el[0]] = sim[el]

        # shift all values so there are no negative numbers
        if np.amin(sim_matrix) < 0:
            sim_matrix += -1*np.amin(sim_matrix)

        clustering = SpectralClustering(n_clusters=n_clusters, affinity="precomputed").fit(sim_matrix)
        labels = clustering.labels_

        return labels

    def calculate_composition(self,collate=False):
        """
        Function that returns the frequency of each letter in the seq_set.
        
        If collate is True, a single frequency dictionary is returned for the entire
        sequence set. If False, a list of frequency dictionaries is returned, corresponding
        to each sequence in the set.
        """
        if not collate:
            all_freqs = []
            for seq in self.records:
                freq = seq.calculate_composition()
                all_freqs.append(freq)

            return all_freqs

        else:
            for i, seq in enumerate(self.records):
                if i == 0:
                    freq = seq.calculate_composition()
                else:
                    freq = add_dicts(freq, seq.calculate_composition())

            return freq

    def filter_by_weight(self,threshold=1000,remove_below=True):
        """Function to filter sequences in seq_set by weight (in units of Da)"""
        new_recs = []
        for seq in self.records:
            weight = seq.calculate_weight()
            if (weight >= threshold and remove_below) or (weight <= threshold and not remove_below):
                new_recs.append(seq)

        self.records = new_recs