"""
This file holds the sequence set (seq_set) class and attributed methods.
Includes I/O of .fasta files, and various filtering methods.
"""

import time
import copy
from sklearn.cluster import SpectralClustering
import numpy as np
import networkx as nx
from Bio import SeqIO
from sequence_analysis.utils import add_dicts
from sequence_analysis.sequence import sequence
from sequence_analysis.utils import dna_alphabet, gen_non_overlapping_points, generate_random_color
from sequence_analysis.utils import rna_alphabet, diff_letters, write_in_columns, aa_alphabet
from sequence_analysis.pairwise_alignment import pairwise_alignment
import re

class seq_set:
    """This class holds a list of sequence objects of a given type."""
    def __init__(self, list_of_sequences=None, file_name=None, seq_type=None):
        if list_of_sequences is None:
            list_of_sequences = []
        self.records = list_of_sequences
        self.type = seq_type

        self.sim_matrix = None
        self.cluster_labels = None

        if len(list_of_sequences) == 0 and file_name is None:
            print("WARNING: sequence set initialized without sequences or file name to read from. Please use read_fasta() to read in sequences or use add_sequence() to add sequences.")

        if file_name is not None:
            self.read_fasta(file_name)

        if seq_type is None and len(self.records) > 0:
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

    def __getitem__(self, key):
        """Function to slice into the seq_set object."""
        if isinstance(key, slice):
            new_set = copy.copy(self)
            new_set.records = new_set.records[key]
            return new_set
        elif isinstance(key, int):
            return copy.copy(self.records[key])
        else:
            raise TypeError

    def get_len(self):
        """Function that returns the number of sequences within seq_set."""
        # TODO: switch to len() in all instances and get rid of this function.
        return len(self.records)

    def write_fasta(self, file_name):
        """"Function to write all sequences within seq_set to a fasta file."""
        file = open(file_name, 'w', encoding="utf-8")
        for i, seq in enumerate(self.records):
            seq_string = seq.seq
            if seq.name is not None:
                name = ">" + seq.name
            else:
                name = f"> {i}"
            file.write(name)
            file.write('\n')
            write_in_columns(file, seq_string, ncols=79)
        file.close()

    def write_fasta_in_parts(self, file_name, n_seq=1000):
        """
        Function to write N fasta files 
        file_name_0.fasta
        file_name_1.fasta
        ...
        file_name_N.fasta

        such that each file contains n_seq sequences.
        """
        # strip any file extention
        # TODO: make this more robust?
        file_name = file_name.split('.')[0]

        n_seqs = len(self)
        if n_seqs % n_seq == 0:
            number_of_files = len(self) // n_seq
        else:
            number_of_files = len(self) // n_seq + 1
        for file_num in range(number_of_files):
            start = file_num * n_seq
            end = (file_num + 1) * n_seq
            curr_file_name = f"{file_name}_{file_num}.fasta"
            curr_set = seq_set(list_of_sequences=self.records[start:end])
            curr_set.write_fasta(curr_file_name)

        return

    def filter_by_six_frame_check_pattern(
            self, regex, overwrite_frame_shifted=True,
            min_orf_len=90):
        """
        Function that checks pattern in the translated sequence and returns
        the translated protein sequence.

        If overwrite_frame_shifted is True, sequences in the set will be
        overwritten by the correct frame-shifted rna or dna sequence.
        """
        assert (self.type in ['dna', 'rna'])

        true_records = []
        untranslated_records = []
        for seq in self.records:
            true, untranslated = seq.six_frame_check(regex, min_orf_len=min_orf_len)
            if true is not None:
                true_records.append(true)
                untranslated_records.append(untranslated)

        if overwrite_frame_shifted:
            self.records = untranslated_records

        return true_records

    def remove_duplicates(self, alphabetize=True):
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
            for i, record in enumerate(self.records):
                if i == 0:
                    unique_records.append(self.records[0])
                else:
                    if record != self.records[i - 1]:
                        unique_records.append(self.records[i])
                    else:
                        unique_records[-1].name += f"_{record.name}"

            self.records = unique_records

    def get_frequencies(self):
        """Function that gets frequencies of all unique sequences."""
        records = sorted(self.records)

        unique_records = []
        frequencies = []
        for i, record in enumerate(records):
            if i == 0:
                count = 1
                unique_records.append(records[0])
            else:
                if record == records[i - 1]:
                    count += 1
                else:
                    frequencies.append(count)
                    count = 1
                    unique_records.append(record)

        frequencies.append(count)
        return unique_records, frequencies

    def sort_by_frequency(self, increasing=False):
        """Function to sort sequences in decreasing frequencies."""
        unique_records, frequencies = self.get_frequencies()

        if increasing:
            order = np.argsort(frequencies)
        else:
            order = np.argsort(frequencies)[::-1]
        sorted_seqs = [unique_records[i] for i in order]
        sorted_freqs = [frequencies[i] for i in order]

        sorted_seq_set = seq_set(list_of_sequences=sorted_seqs)

        return sorted_seq_set, sorted_freqs

    def alphabetize(self):
        """Function to sort sequences within seq_set in alphabetical order."""
        records = sorted(self.records)
        self.records = records

    def get_letters(self):
        """"Function to get all unique letters in the seq_set sequences."""
        all_letters = set()
        for seq in self.records:
            seq_string = seq.seq.upper()
            all_letters = all_letters.union(set([*seq_string]))
        return all_letters

    def set_type(self, seq_type=None):
        """Function to detect and set the type of sequences in seq_set."""
        # TODO: do further testing of type detection
        all_letters = self.get_letters()
        if 'N' in all_letters:
            all_letters.remove('N')

        if seq_type is not None:
            self.type = seq_type
            # TODO: check that the assigned type is consistent
            return

        if len(diff_letters.intersection(all_letters)):
            self.type = 'protein'
        elif set(dna_alphabet).issubset(all_letters) or set(rna_alphabet).issubset(all_letters):
            if 'U' in all_letters:
                self.type = 'rna'
            else:
                # TODO: this could fail if we have RNA sequences that happens
                # to have no Us, but that's unlikely
                self.type = 'dna'
        else:
            print("I can't assign a type, please specify manually.")

        # translate the type to all sequences
        for i in range(len(self.records)):
            self.records[i].type = self.type

    def read_fasta(self, file_name):
        """Function to read sequences into seq_set from a .fasta file."""
        if len(self.records) > 0:
            print("Warning: overwriting existing data")
            self.records = []
        for record in SeqIO.parse(file_name, "fasta"):
            seq = sequence(str(record.seq), record.name)
            self.records.append(seq)

        self.set_type()

        for seq in self.records:
            seq.type = self.type

    def add_sequence(self, sequences):
        """Function that adds another sequence to seq_set."""
        # TODO: refactor, because you have an add_set method as well
        if isinstance(sequences, list):
            for seq in sequences:
                assert (isinstance(seq, sequence))
                self.records.append(seq)
            old_type = self.type
            self.set_type()
            if old_type is not None:
                assert (old_type == self.type)
            self.records[-1].type = self.type
        elif isinstance(sequences, sequence):
            self.records.append(sequences)
            old_type = self.type
            self.set_type()
            if old_type is not None:
                assert (old_type == self.type)
            self.records[-1].type = self.type
        else:
            print("ERROR: incorrect format for adding sequence.",
                  "Please provide a list of sequences or a sequence object. Sequence was not added.")

    def add_set(self, set2):
        """Function that merges two sets of sequences (two seq_sets)."""
        # make sure the set being added is of type seq_set
        if not isinstance(set2, seq_set):
            print("ERROR: cannot add the two sets. Set2 is not of type seq_set.")
            return

        # make sure sequence types are the same
        if set2.type != self.type:
            print("ERROR: cannot add two sets. They are not of the same type.")
            return

        # add set2's records to the list of records
        self.records += set2.records

    def filter_by_pattern(self, regex):
        """Filter sequences in seq_set by the given regex pattern."""
        new_records = []
        for seq in self.records:
            if seq.check_for_pattern(regex):
                new_records.append(seq)

        self.records = new_records

    def get_similarity_matrix(self,
                              algorithm="biopython-global",
                              use_blosum_50=False,
                              match=2,
                              unmatch=-1,
                              gap=-0.1,
                              gap_open=-0.5,
                              verbose=False):
        """
        Function that calculates similarity matrix for the sequences in seq_set.

        For a seq_set containing N sequences, calculates (N-1)N/2 pairwise alignment scores.
        """

        # pairwise alignment between all pairs of sequences
        n_sequences = self.get_len()
        similarity_matrix = {}
        for i in range(n_sequences):
            for j in range(i + 1, n_sequences):
                start = time.time()
                alignment = pairwise_alignment(self.records[i],
                                               self.records[j],
                                               algorithm=algorithm,
                                               use_blosum_50=use_blosum_50,
                                               match=match,
                                               unmatch=unmatch,
                                               gap=gap,
                                               gap_open=gap_open)
                alignment.align()
                similarity_matrix[(i, j)] = alignment.score
                end = time.time()

                if verbose:
                    print(i, j, end - start)

        self.sim_matrix = similarity_matrix

    def filter_by_frequency(self, threshold=10):
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
            seq.remove_before_pattern(regex, verbose=verbose)

    def cluster(self,
                n_clusters,
                algorithm="biopython-global",
                use_blosum_50=False,
                match=2,
                unmatch=-1,
                gap=-0.1,
                gap_open=-0.5,
                verbose=False):
        """
        Function that clusters the sequences in seq_set.

        Similarity matrix is calculated based on pairwise alignments, and then used to
        perform spectral clustering.
        """
        if self.sim_matrix is None:
            self.get_similarity_matrix(algorithm=algorithm,
                                       use_blosum_50=use_blosum_50,
                                       match=match,
                                       unmatch=unmatch,
                                       gap=gap,
                                       gap_open=gap_open,
                                       verbose=verbose)
        else:
            n_sequences = len(self)
            n_sim_matrix_el = len(self.sim_matrix)
            if n_sim_matrix_el != (n_sequences - 1) * n_sequences / 2:
                self.get_similarity_matrix(algorithm=algorithm,
                                           use_blosum_50=use_blosum_50,
                                           match=match,
                                           unmatch=unmatch,
                                           gap=gap,
                                           gap_open=gap_open,
                                           verbose=verbose)

        sim = self.sim_matrix

        n_sequences = len(self)
        sim_matrix = np.zeros((n_sequences, n_sequences))
        for el in sim.keys():
            sim_matrix[el[0], el[1]] = sim[el]
            sim_matrix[el[1], el[0]] = sim[el]

        # shift all values so there are no negative numbers
        if np.amin(sim_matrix) < 0:
            sim_matrix += -1 * np.amin(sim_matrix)

        clustering = SpectralClustering(
            n_clusters=n_clusters,
            affinity="precomputed").fit(sim_matrix)
        self.cluster_labels = clustering.labels_

    def visualize_clusters(self):
        """Function that visualizes clusters."""
        # TODO: refactor (break into more functions? too many local vars)
        if self.cluster_labels is None:
            print(
                "ERROR: cannot visualize clusters because they don't exist. Run cluster() first.")
            return

        assert (self.sim_matrix is not None)

        n_sequence = len(self)

        n_clusters = np.unique(self.cluster_labels).size

        G = nx.complete_graph(n_sequence)
        colors = generate_random_color(n_clusters)
        node_colors = []
        subgraph_indices = {i: [] for i in range(n_clusters)}
        for i in range(n_sequence):
            label = self.cluster_labels[i]
            subgraph_indices[label].append(i)
            node_colors.append(colors[label])

        # generate positions for nodes
        pos = {}
        box_lengths = []
        for i in range(n_clusters):
            subgraph = G.subgraph(subgraph_indices[i])
            pos_subgraph = nx.spring_layout(subgraph)

            coords = np.array([i for _, i in pos_subgraph.items()])

            center = np.mean(coords, axis=0)

            span = [[np.amin(coords[:, 0]), np.amax(coords[:, 0])],
                    [np.amin(coords[:, 1]), np.amax(coords[:, 1])]]

            box_lengths.append(
                max(abs(span[0][1] - span[0][0]), abs(span[1][1] - span[1][0])))

            for _, coord in pos_subgraph.items():
                coord[0] -= center[0]
                coord[1] -= center[1]

            pos |= pos_subgraph

        # translate positions of each cluster
        separation = max(box_lengths)
        radius = 4 * separation / (2 * np.pi) * n_clusters
        cluster_centroids = []
        for cluster in range(n_clusters):
            angle = cluster * (2 * np.pi) / n_clusters
            centroid_x = radius * np.cos(angle)
            centroid_y = radius * np.sin(angle)
            cluster_centroids.append([centroid_x, centroid_y])

        for ind, coord in pos.items():
            label = self.cluster_labels[ind]
            coord[0] += cluster_centroids[label][0]
            coord[1] += cluster_centroids[label][1]

        nx.draw_networkx_nodes(G, pos=pos, node_color=node_colors)
        # TODO: implement edge weights
        nx.draw_networkx_edges(G, pos=pos)

    def calculate_composition(self, collate=False):
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

    def filter_by_weight(self, threshold=1000, remove_below=True):
        """Function to filter sequences in seq_set by weight (in units of Da)"""
        # TODO: clean up
        new_recs = []

        if remove_below:
            threshold *= -1
        for seq in self.records:
            weight = seq.calculate_weight() * -1
            # if (weight >= threshold and remove_below) or (weight <= threshold
            # and not remove_below):
            if weight <= threshold:
                new_recs.append(seq)

        self.records = new_recs

    def crop_first_N(self, N):
        """Function to cut sequences such that the first N chars are preserved."""
        for seq in self.records:
            n_cut = min(len(seq), N)
            seq.seq = seq.seq[:n_cut]

    def add_gaps_to_align_pattern(self, pattern):
        """
        Given a pattern, add gaps to the beginning of the sequence so that the patterns
        align for all the sequences.

        Takes into account only the first occurence of the pattern.
        """

        pattern = re.compile(pattern)
        indices = []
        for seq in self:
            try:
                idx = pattern.search(seq.seq).start()
            except AttributeError:
                idx = -1
            indices.append(idx)

        max_idx = max(indices)
        for i, seq in enumerate(self.records):
            idx = indices[i]
            if idx != -1:
                n_chars = max_idx - idx
                seq.seq = '-' * n_chars + seq.seq

        return

    def create_sequence_logo(self, char_range):
        """
        Function to generate the information necessary for a sequence logo.
        """
        if self.type == 'dna':
            alphabet = dna_alphabet
            s = 4
        elif self.type == 'rna':
            alphabet = rna_alphabet
            s = 4
        elif self.type == 'protein':
            alphabet = aa_alphabet
            s = 20
        else:
            print("ERROR: sequence set type not detected.")
            raise TypeError

        # create frequency matrix
        start = char_range[0]
        end = char_range[1]
        frequencies = [[0 for _ in range(start, end)] for _ in alphabet]

        # create matrix that will hold the number of characters at that position
        counts = [0 for _ in range(start, end)]

        for seq in self:
            for i, char in enumerate(seq):
                if i >= start and i < end:
                    idx = i - start
                    counts[idx] += 1
                    char_idx = alphabet.index(char)
                    frequencies[char_idx][idx] += 1

        # normalize freqs
        for l, ff in enumerate(frequencies):
            for i, _ in enumerate(ff):
                frequencies[l][i] /= counts[i]

        # TODO: cleanup
        # calculating the information content
        e_n = [1./np.log(2) * (s-1.)/(2*c) for c in counts]
        R_i = [0 for _ in range(start, end)]
        H_i = [0 for _ in range(start, end)]
        height = [[0 for _ in range(start, end)] for _ in alphabet]
        for i in range(len(alphabet)):
            for j in range(len(counts)):
                freq = frequencies[i][j]
                if freq != 0:
                    H_i[j] += -1 * freq * np.log2(freq)

        for i in range(len(counts)):
            R_i[i] = np.log2(s) - (H_i[i] + e_n[i])

        for i in range(len(alphabet)):
            for j in range(len(counts)):
                height[i][j] = frequencies[i][j] * R_i[j]

        return alphabet, height