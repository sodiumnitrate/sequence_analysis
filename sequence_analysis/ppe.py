"""
Contains the PPE class that generates an alphabet of frequent words.
Can be further used for bag-of-words type of applications.

Algorithm described in Asgari et al., Scientific Reports, 2019.

Algorithmic complexity: O(Nm^2) with N=number of seqs, m=max seq. length,
if max_iter was inf.
(Happens in the case where sequences are completely dissimilar and have
absolutely no repetitions within themselves either. A more reasonable 
expectation is something like O(Nmlog(m)).)

TODO:
- currently duplicates the entire seq_set data. consider keeping things
as sequences and doing operations on the fly.
- should be able to simplify alphabet update.
- check for list of list of strings with single strings?
"""

import copy
import numpy as np

from sklearn.cluster import SpectralClustering

from sequence_analysis.seq_set import seq_set
from sequence_analysis.utils import find_kmers_in_list, merge_dicts
from sequence_analysis.sequence import sequence


class PPE:
    def __init__(self, sequences):
        """
        Function to initialize the PPE.
        """
        self.sequences = sequences

        self.n_sequences = None
        self.alphabet = None
        self.merge_operations = []
        self.frequency_of_seqs_with_most_common = 0
        self.frequency_of_most_common_motif = 0
        self.n_iter = 0
        self.motifs = None
        self.converged = None
        self.sim_matrix = None

        self.process_sequences()
        self.generate_alphabet()

        self.vectors = None
        self.adjacency = None

    def process_sequences(self):
        """
        Convert sequences to list of characters.
        
        Accepted input types:
        - seq_set object
        - list of strings
        - list of lists that contain chars
        """
        if isinstance(self.sequences, seq_set):
            new_sequences = []
            for seq in self.sequences:
                new_sequences.append(list(seq.seq))
            self.sequences = new_sequences
        elif isinstance(self.sequences, list):
            # sequences is a list
            list_of_strings = all([isinstance(seq, str) for seq in self.sequences])
            list_of_list_of_strings = all([all([isinstance(s, str) for s in seq]) for seq in self.sequences])
            list_of_list_of_strings = list_of_list_of_strings and all([isinstance(seq, list) for seq in self.sequences])
            list_of_sequences = all([isinstance(seq, sequence) for seq in self.sequences])
            if list_of_strings:
                new_sequences = [list(seq) for seq in self.sequences]
                self.sequences = new_sequences
            elif list_of_sequences:
                # list of sequence objects
                new_sequences = [list(seq.seq) for seq in self.sequences]
                self.sequences = new_sequences
            elif list_of_list_of_strings:
                # already in the form we want
                pass
            else:
                print(f"ERROR: sequence list of {type(self.sequences[0])} is not supported.")
                raise TypeError
        else:
            print(f"ERROR: sequences of type {type(self.sequences)} not supported.")
            raise TypeError

        self.n_sequences = len(self.sequences)

    def generate_alphabet(self):
        """
        Given a list of sequences, generate an alphabet.
        """
        alphabet = set()
        for seq in self.sequences:
            alphabet = alphabet.union(set(seq))

        self.alphabet = alphabet

    def merge_and_update_alphabet(self):
        """
        Do a pass at merging twomers.
        """
        self.frequency_of_seqs_with_most_common = 0
        new_sequences = []
        for seq in self.sequences:
            symbol_in_seq = False
            new_seq = []
            n = len(seq)
            i = 0
            while i < n - 1:
                curr_tuple = (seq[i], seq[i+1])
                if curr_tuple == self.merge_operations[-1]:
                    new_seq.append(seq[i] + seq[i+1])
                    i += 1
                    symbol_in_seq = True
                else:
                    new_seq.append(seq[i])

                i += 1

            if i == n - 1:
                new_seq.append(seq[-1])

            new_sequences.append(copy.copy(new_seq))
            if symbol_in_seq:
                self.frequency_of_seqs_with_most_common += 1

        # TODO: should be able to do this on the fly as well. not sure worth it.
        self.sequences = new_sequences
        self.generate_alphabet()

    def find_most_common_twomer(self):
        """
        Function that finds the most common twomer in the list of sequences.
        """
        all_twomers = {}
        for seq in self.sequences:
            curr_twomers = find_kmers_in_list(seq, 2)
            all_twomers = merge_dicts(all_twomers, curr_twomers)
        most_freq = max(all_twomers.items(), key=lambda x:x[1])
        self.merge_operations.append(most_freq[0])
        self.frequency_of_most_common_motif = most_freq[1]

    def iterate(self, max_iter=100000):
        """
        Function to iteratively find motifs.
        Stopping criterion:
        self.frequency_of_seqs_with_most_common >= self.frequency_of_most_common_motif

        or self.n_iter = max_iter
        """

        self.frequency_of_seqs_with_most_common = 0
        self.frequency_of_most_common_motif = 1
        self.n_iter = 0
        while self.frequency_of_seqs_with_most_common < self.frequency_of_most_common_motif and self.n_iter < max_iter:
            self.find_most_common_twomer()
            self.merge_and_update_alphabet()
            self.n_iter += 1

        if self.frequency_of_seqs_with_most_common < self.frequency_of_most_common_motif:
            self.converged = False
        else:
            self.converged = True

        self.generate_list_of_found_motifs()

    def generate_list_of_found_motifs(self):
        """
        Function to generate a list of found motifs from the current alphabet.
        """
        alpha_counts = {a: 0 for a in self.alphabet}
        for seq in self.sequences:
            for a in seq:
                alpha_counts[a] += 1
        self.motifs = sorted(alpha_counts.items(), key=lambda x:x[1], reverse=True)

    def generate_vectors(self):
        """
        Once the motifs have been determined, assing a vector to each sequence that
        corresponds to the number of occurences of each letter in the alphabet.
        """
        if self.converged is None:
            print("ERROR: run iterate() first.")
            return

        self.vectors = [[0 for _ in self.alphabet] for _ in self.sequences]
        for j, seq in enumerate(self.sequences):
            for i, a in enumerate(self.alphabet):
                ct = seq.count(a)
                self.vectors[j][i] = ct

    def find_distance_between_seqs(self, i, j):
        """
        Once the motifs are found and bag-of-words vectors are calculated,
        calculate the Euclidean distance between two sequences.
        """
        if self.vectors is None:
            print("ERROR: vectors are not calculated!")
            return

        vec1 = np.array(self.vectors[i])
        vec2 = np.array(self.vectors[j])

        # TODO: operate by dist^2 instead?
        return np.linalg.norm(vec1-vec2)

    def adjacency_list(self):
        """
        Once the motifs are found and vectors calculated, get an adjacency
        list between all pairs of sequences.
        """
        self.adjacency = {}

        for i in range(self.n_sequences):
            for j in range(i+1, self.n_sequences):
                self.adjacency[(i,j)] = self.find_distance_between_seqs(i,j)

    def cluster(self, n_clusters):
        """
        Function to cluster using spectral clustering, using the adjacency list.
        """
        sim = self.adjacency
        sim_matrix = np.zeros((self.n_sequences, self.n_sequences))
        for el, val in sim.items():
            sim_matrix[el[0], el[1]] = val
            sim_matrix[el[1], el[0]] = val

        # sim_matrix contains distances. Normalize by max dist, and subtract from 1 to get similarity
        sim_matrix = np.array(sim_matrix) / np.amax(sim_matrix)
        sim_matrix = 1 - sim_matrix

        self.sim_matrix = sim_matrix

        clustering = SpectralClustering(n_clusters=n_clusters,
                                        affinity="precomputed").fit(sim_matrix)

        return clustering