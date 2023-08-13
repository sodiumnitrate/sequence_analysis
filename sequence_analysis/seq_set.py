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
from sequence_analysis.utils import merge_dicts

class seq_set:
    """This class holds a list of sequence objects of a given type."""
    def __init__(self, list_of_sequences=None, file_name=None, seq_type=None):
        if list_of_sequences is None:
            list_of_sequences = []
        self.records = list_of_sequences
        self.type = seq_type

        self.sim_matrix = None
        self.cluster_labels = None

        self.name_dict = None

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

    def __add__(self, other):
        """Overwrites __add__ so that two seq_sets can be merged."""
        if not isinstance(other, seq_set):
            print("ERROR: both inputs must be of type seq_set.")
            raise TypeError

        if self.type != other.type:
            print(f"ERROR: can't add a set with type {self.type} and {other.type}.")
            raise TypeError

        new_set = seq_set(list_of_sequences=self.records+other.records)
        new_set.type = self.type
        return new_set

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
                        if record.name not in unique_records[-1].name:
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

    def filter_by_pattern(self, regex):
        """Filter sequences in seq_set by the given regex pattern."""
        new_records = []
        for seq in self.records:
            if seq.check_for_pattern(regex):
                new_records.append(seq)

        self.records = new_records

    def unalign(self):
        """
        Removes gaps from all sequences so that the sequence set is unaligned.
        """
        new_recs = []
        for rec in self.records:
            seq_str = rec.seq.replace('-','')
            seq = sequence(seq_str)
            seq.name = rec.name
            seq.type = rec.type
            new_recs.append(seq)

        self.records = new_recs

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
                if char == '-':
                    continue
                if i >= start and i < end:
                    idx = i - start
                    counts[idx] += 1
                    char_idx = alphabet.index(char)
                    frequencies[char_idx][idx] += 1

        # normalize freqs
        for l, ff in enumerate(frequencies):
            for i, _ in enumerate(ff):
                if counts[i] != 0:
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

    def find_kmers(self, k):
        """
        Function to find kmers of length k in all the sequences in the set, with
        their frequencies.
        """

        kmers = {}
        for seq in self:
            curr_kmer = seq.find_kmers(k)
            kmers = merge_dicts(kmers, curr_kmer)

        return kmers
    
    def find_subset_with_names(self, names, suppress_warning=False):
        """
        Given a list of names, return a seq_set with sequences that have
        matching names.
        """
        # make sure we have a unique list of names
        names = list(set(names))

        # get names of seqs in the sequence set
        if self.name_dict is None:
            # makes subsequent searches faster
            self_names = [s.name for s in self]
            if len(list(set(self_names))) != len(self_names):
                print("WARNING: there are multiple sequences with the same name. I will take the last one.")

            self.name_dict = {s.name:i for i, s in enumerate(self)}

        subset = []
        for name in names:
            try:
                idx = self.name_dict[name]
            except KeyError:
                if not suppress_warning:
                    print(f"WARNING: no sequence with name {name} found.")
                continue
            seq = self.records[idx]
            subset.append(seq)

        if len(subset) == 0:
            return None

        return seq_set(subset)