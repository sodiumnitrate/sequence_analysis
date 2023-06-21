"""
This file holds the GibbsSampler class, that, given a seq_set,
sets up a Gibbs sampler to find motifs.

Implemented from Lawrence et al. Science 1993
"""
import random
import copy
import numpy as np
from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence
from sequence_analysis.utils import aa_alphabet, dna_alphabet, rna_alphabet

import pdb

class GibbsSampler:
    def __init__(self, sequences, word_length):
        """
        Function to initialize GibbsSampler object.
        """
        self.sequences = sequences
        self.word_length = word_length

        self.alphabet = None
        self.alphabet_dict = None
        self.pseudocounts = None

        self.check_input_and_set_alphabet()

        self.indices = None
        self.pssm = None
        self.prev_indices = None
        self.background_probabilities = None

        self.total_score = None
        self.converged = None
        self.n_iter = None

    def check_input_and_set_alphabet(self):
        """
        Function to check input data and choose the appropriate alphabet.
        """
        if not isinstance(self.word_length, int):
            print("ERROR: word_length must be an integer.")
            raise TypeError
        if not isinstance(self.sequences, seq_set):
            if isinstance(self.sequences, list):
                if isinstance(list[0], str):
                    # we have a list of strings, convert to seq
                    sequences = [sequence(s) for s in self.sequences]
                    self.sequences = seq_set(list_of_sequences=sequences)
                else:
                    print("ERROR: sequences must be a list of strings or a seq_set object.")
                    raise TypeError
            else:
                print("ERROR: sequences must be a list of strings or a seq_set object.")
                raise TypeError

        # check that word_length is shorter than the shortest sequence
        min_length = min([len(seq) for seq in self.sequences])
        if self.word_length > min_length:
            print(f"ERROR: word length, {self.word_length}, can't be larger than the minimum sequence length, {min_length}.")

        # check seq types
        if self.sequences.type == 'protein':
            self.alphabet = aa_alphabet
        elif self.sequences.type == 'rna':
            self.alphabet = rna_alphabet
        elif self.sequences.type == 'dna':
            self.alphabet = dna_alphabet
        else:
            print(f"ERROR: sequences of type {self.sequences.type} cannot be analyzed.")
            raise TypeError

        # set alphabet_dict
        self.alphabet_dict = {char:i for i, char in enumerate(self.alphabet)}
        self.set_pseudocounts()

    def set_pseudocounts(self):
        """
        Set pseudocounts (see Lawrence et al., Science 1993).

        Currently just ones.
        """
        self.pseudocounts = [1 for _ in self.alphabet]

    def initialize_sampler(self, random_init=True, indices=None):
        """
        Function to initialize starting point.
        """
        if random_init:
            self.indices = [random.randint(0, len(seq) - self.word_length) for seq in self.sequences]
        else:
            if not indices:
                print("ERROR: random_init is set to False but no indices provided.")
                return
            else:
                if not isinstance(indices, list) or not isinstance(indices[0], int):
                    print("ERROR: start indices must be a 1d list of integers.")
                    return
                if len(indices) != len(self.sequences):
                    print(f"ERROR: list of start indices have length {len(indices)} but we have {len(self.sequences)} sequences.")
                    return
                self.indices = indices

    def init_pssm(self):
        """
        Function to initialize PSSM (to zeros).
        """
        if self.alphabet is None:
            print("ERROR: can't initialize PSSM -- alphabet is not set!")
            return

        self.pssm = [[0 for _ in range(self.word_length)] for _ in self.alphabet]

    def calculate_background_probabilities(self, seq_idx):
        """
        Calculate background probabilities of letters outside the current
        chosen regions.

        seq_idx is a list of sequence indices.
        """
        # TODO: make this position-specific
        if self.alphabet is None:
            print("ERROR: can't calculate background probabilities -- alphabet is not set.")
            return
        
        background_probabilities = [0 for _ in self.alphabet]
        freq = 0
        for idx in seq_idx:
            start = self.indices[idx]
            end = start + self.word_length
            for i, char in enumerate(self.sequences[idx]):
                if i < start or i >= end:
                    background_probabilities[self.alphabet_dict[char]] += 1
                    freq += 1

        B = sum(self.pseudocounts)
        self.background_probabilities = [(p + self.pseudocounts[i])/(freq + B) for i, p in enumerate(background_probabilities)]

    def calculate_pssm(self, exclude=None):
        """
        Calculate PSSM for the current set of indices.
        """
        seq_idx = list(range(len(self.sequences)))
        if exclude:
            if not isinstance(exclude, int):
                print("ERROR: exclude must be an integer, if not None.")
                raise TypeError
            seq_idx.remove(exclude)

        self.init_pssm()
        self.calculate_background_probabilities(seq_idx)
        freqs = [0 for _ in range(self.word_length)]

        for idx in seq_idx:
            start = self.indices[idx]
            end = start + self.word_length
            seq_str = self.sequences[idx].seq
            for i in range(start, end):
                char = seq_str[i]
                self.pssm[self.alphabet_dict[char]][i-start] += 1
                freqs[i-start] += 1

        B = sum(self.pseudocounts)
        for aa_idx in range(len(self.pssm)):
            for pos_idx in range(self.word_length):
                f = freqs[pos_idx]
                self.pssm[aa_idx][pos_idx] += self.pseudocounts[aa_idx]
                self.pssm[aa_idx][pos_idx] /= (f + B)
                bg = self.background_probabilities[aa_idx]
                self.pssm[aa_idx][pos_idx] /= bg
                
    def calculate_profile_score(self, i, j):
        """
        Calculates profile score given the PSSM matrix and background probabilities,
        with sequence i excluded, and start index for i is j.
        """
        # TODO: make sure pssm calculated etc.
        # TODO: how to check if pssm is up-to-date?
        seq = self.sequences[i]
        score = 0
        for idx, char in enumerate(seq[j:j+self.word_length]):
            aa_idx = self.alphabet_dict[char]
            score += np.log(self.pssm[aa_idx][idx])

        return score

    def calculate_overall_score(self):
        """
        Function to calculate overall score. Assumes iteration is run and converged.
        """
        total_score = 0
        for i in range(len(self.sequences)):
            j = self.indices[i]
            total_score += self.calculate_profile_score(i, j)

        self.total_score = total_score

    def iterate(self, max_iter=100):
        """
        Iterate until the indices don't change, or max_iter is reached.
        """
        if not self.indices:
            print("WARNING: start indices not initialized. Setting randomly.")
            self.initialize_sampler()

        self.prev_indices = None
        n_iter = 0
        reached_max_iter = False
        while self.prev_indices != self.indices:
            self.prev_indices = copy.copy(self.indices)

            # iterate through every string
            for i in range(len(self.sequences)):
                # calculate pssm for all seqs except i
                self.calculate_pssm(exclude=i)

                # find index in i where the profile matches best
                best_score = None
                for j in range(len(self.sequences[i]) - self.word_length + 1):
                    score = self.calculate_profile_score(i, j)
                    # TODO: implemented random weighted selection, not steepest descent
                    if best_score is None:
                        best_score = score
                        best_pos = j
                    elif score > best_score:
                        best_score = score
                        best_pos = j

                self.indices[i] = best_pos
            n_iter += 1
            if n_iter >= max_iter:
                print(f"WARNING: stopping because max_iter={max_iter} is reached.")
                reached_max_iter = True

        if not reached_max_iter:
            print(f"GibbsSampler iterator converged in {n_iter} steps.")
            self.converged = True
        else:
            self.converged = False

        self.n_iter = n_iter
        self.calculate_pssm()
        self.calculate_overall_score()
        self.shift_all_indices()
        print(f"Final score is: {self.total_score}.")

    def shift_all_indices(self):
        """
        Function to shift all indices and recalculate scores to solve the "phase problem."
        (See Lawrence et al., Science 1993)
        """
        best_indices = copy.copy(self.indices)
        best_score = self.total_score
        original_indices = copy.copy(self.indices)
        for i in range(-1*self.word_length, self.word_length + 1):
            new_indices = copy.copy(original_indices)
            for j in range(len(new_indices)):
                new_indices[j] += i

            # check indices
            valid_idx = True
            for j, idx in enumerate(new_indices):
                if idx < 0:
                    valid_idx = False
                    break
                if idx + self.word_length >= len(self.sequences[j]):
                    valid_idx = False
                    break

            if not valid_idx:
                continue

            self.indices = new_indices
            self.calculate_pssm()
            self.calculate_overall_score()
            if self.total_score > best_score:
                best_score = self.total_score
                best_indices = copy.copy(self.indices)

        self.indices = best_indices
        self.calculate_pssm()
        self.calculate_overall_score()
        assert best_score == self.total_score

    def print_motifs(self):
        """
        Print the detected motifs.
        """
        if self.indices is None:
            print("ERROR: did you initialize and run .iterate()?")
            return
        
        for i in range(len(self.indices)):
            seq = self.sequences[i]
            idx = self.indices[i]
            print(seq[idx:idx+self.word_length])
