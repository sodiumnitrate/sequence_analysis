"""
This file holds the sequence class and related methods.
"""
from sequence_module import PySequence

class Sequence:
    """This class corresponds to a single biological sequence (protein, rna, or dna)"""
    def __init__(self, seq, name=None, seq_type=None):
        """This function initializes the sequence object."""
        self.py_seq = PySequence(seq, name)
        if seq_type is None:
            self.py_seq.set_type("")
        else:
            self.py_seq.set_type(seq_type)

    def __eq__(self, other):
        """
        This function overloads the __eq__ function such that if two sequences are of the
        same type and have the same exact string for sequence, they are deemed equal.
        """
        if isinstance(self, type(other)):
            if self.py_seq.type == other.py_seq.type:
                return self.py_seq.get_sequence() == other.py_seq.get_sequence()
            else:
                return False
        else:
            return False

    def __str__(self):
        """Function that overloads __str__ for sequence object."""
        return f"{self.py_seq.get_sequence()}"

    def __repr__(self):
        """Function that overloas __repr__ for sequence object."""
        return f"<Sequence object of type {self.py_seq.get_type()} with '{self.py_seq.get_sequence():.5}...' at {hex(id(self))}>"

    def __len__(self):
        """Overwrites __len__ to return the number of chars in seq.seq"""
        return self.py_seq.length()

    def __getitem__(self, key):
        """Overwrites __getitem__ so that a part of the sequence can be returned."""
        if isinstance(key, slice) or isinstance(key, int):
            return self.py_seq.get_sequence()[key]
        else:
            raise TypeError

    def set_type(self, seq_type=None):
        """Function that sets the type of sequence based on the letters it contains."""
        if seq_type is None:
            self.py_seq.set_type("")
        else:
            self.py_seq.set_type(seq_type)

    def translate(self):
        """
        Function that translates DNA or RNA sequence into protein.
        """
        if self.py_seq.get_type() != 'rna' and self.py_seq.get_type() != 'dna':
            print(f"WARNING: no translation for sequence type {self.py_seq.get_type()}.")
            return
        '''
        s = Seq(self.seq.replace('-',""))
        s = str(s.translate())
        s = s.strip("X")

        return sequence(s,name=self.name,seq_type='protein')
        '''
        s_cpp = self.py_seq.translate()
        s = Sequence(s_cpp.get_sequence(), s_cpp.get_name(), s_cpp.get_type())
        return s

    def find_codon(self, codon):
        """
        Function to search for codons. Returns index of first match.
        Returns None if not rna or dna.
        Returns -1 if no match.
        """
        return self.py_seq.find_codon(codon)

    def reverse_complement(self):
        if self.py_seq.get_type() != "rna" and self.py_seq.get_type() != "dna":
            print(f"ERROR: cannot reverse complement: the sequence is of type {self.py_seq.get_type()}.")
            return None

        seq = self.seq

        seq = seq.replace("C","X")
        seq = seq.replace("G","C")
        seq = seq.replace("X","G")

        if self.type == "rna":
            seq = seq.replace("A","X")
            seq = seq.replace("U","A")
            seq = seq.replace("X","U")
        else:
            seq = seq.replace("A","X")
            seq = seq.replace("T","A")
            seq = seq.replace("X","T")

        # reverse
        seq = seq[::-1]
        return sequence(seq, seq_type=self.type, name=self.name)

    def frame_shift(self, frame=0):
        """Function to frame shift a dna or rna sequence."""
        if self.type != 'rna' and self.type != 'dna':
            print(
                f"ERROR: type {self.type} is not meaningful for frame shifting. Please provide RNA or DNA sequences.")
            return None

        seq = self.seq

        # returns frame shifted sequence
        if frame == 0:
            if len(seq) % 3 == 0:
                seq = sequence(seq)
                seq.type = self.type
                return seq
        elif frame == 1:
            seq = seq[1:]
        elif frame == 2:
            seq = seq[2:]
        else:
            print("ERROR: problem with frame shifting. Frame should be 0, 1, or 2.")
            return None

        if len(seq) % 3 == 0:
            N_to_be_added = 0
        else:
            N_to_be_added = 3 - len(seq) % 3
        seq = seq + "N" * N_to_be_added

        seq = sequence(seq)
        seq.type = self.type

        return seq

    # TODO: DEPRECATE
    def six_frame_check(self, regex,
                              chop_before_first_M=True,
                              min_orf_len=90):
        """
        docstring
        """
        orfs = self.find_open_reading_frames(chop_before_first_M=False, translate=False, min_orf_len=min_orf_len)
        
        if not orfs:
            print("WARNING: no open reading frames were found")

        for orf in orfs:
            orf_seq = sequence(orf, seq_type=self.type)
            seq = orf_seq.translate()
            selected = seq.check_for_pattern(regex)
            if selected:
                ind = max(orf_seq.find_codon(start_codon[0]),
                          orf_seq.find_codon(start_codon[1]))
                if chop_before_first_M:
                    curr = orf_seq.seq[ind:]
                    if len(curr) >= min_orf_len:
                        return sequence(curr, seq_type=self.type).translate().seq, curr
        return None, None

    def find_open_reading_frames(self, chop_before_first_M=True,
                                       translate=True,
                                       min_orf_len=90):
        """
        This function finds all open reading frames. Applies the steps in
        https://vlab.amrita.edu/?sub=3&brch=273&sim=1432&cnt=1
        """
        if self.type != 'rna' and self.type != 'dna':
            print("ERROR: open reading frames can be found for DNA or RNA sequences only.")
            return None
        orfs = []
        for frame in range(3):
            for order in range(2):
                if order == 0:
                    shifted = self.frame_shift(frame=frame)
                else:
                    shifted = self.reverse_complement().frame_shift(frame=frame)

                seq_str = shifted.seq.strip('N')
                n = len(seq_str) // 3
                ptr = 0
                curr_fragment = ""
                for ptr in range(0, 3*n, 3):
                    curr_codon = seq_str[ptr:ptr+3]
                    if curr_codon in stop_codon:
                        if len(curr_fragment) >= min_orf_len:
                            # check if it contains a start codon
                            ind1 = sequence(curr_fragment, seq_type=self.type).find_codon(start_codon[0])
                            ind2 = sequence(curr_fragment, seq_type=self.type).find_codon(start_codon[1])

                            if ind1 != -1 or ind2 != -1:
                                if chop_before_first_M:
                                    curr_fragment = curr_fragment[max(ind1,ind2):]
                                    if len(curr_fragment) < min_orf_len:
                                        continue
                                orfs.append(curr_fragment)
                                curr_fragment = ""
                    else:
                        curr_fragment += curr_codon

        if not translate:
            return orfs

        new_orfs = []
        for orf in orfs:
            seq = sequence(orf, seq_type=self.type).translate()
            new_orfs.append(seq.seq)

        return new_orfs

    # TODO: DEPRECATE
    def check_for_pattern(self, regex):
        """Function that checks for a given regex pattern in the sequence."""
        # TODO: check if regex and seq from the same alphabet?
        return utils.check_for_pattern(self.seq, regex)

    # TODO: DEPRECATE
    def choose_all_matching_patterns(
            self, regex, return_between_matching=False):
        """Function that chooses all matching regex patterns in the sequence."""
        s = self.seq
        p = re.compile(regex)

        x = p.findall(s)

        spans = utils.find_spans(s, x)

        if return_between_matching:
            inverted = utils.invert_spans(spans)
            return x, spans, inverted

        return x, spans

    # TODO: DEPRECATE
    def choose_all_matching_patterns_and_print(self, regex):
        """Function to print spans of matches in blue, and inverted spans in red."""
        matching, spans, inverted = self.choose_all_matching_patterns(
            regex, return_between_matching=True)

        seq = self.seq

        to_print = f""
        initial = seq[:spans[0][0]]
        to_print += initial
        for i, span in enumerate(spans):
            to_print += f"{Fore.BLUE}" + \
                seq[span[0]:span[1] + 1] + f"{Style.RESET_ALL}"
            if i < len(spans) - 1:
                to_print += f"{Fore.RED}" + seq[inverted[i][0]:inverted[i][1] + 1] + f"{Style.RESET_ALL}"
        tail = seq[spans[-1][1] + 1:]
        to_print += tail

        print(to_print)

    # TODO: DEPRECATE
    def remove_before_pattern(self, regex, verbose=True):
        """Function that removes the letters before the matching regex."""
        x, spans = self.choose_all_matching_patterns(regex)
        if verbose:
            if len(spans) == 0:
                print("ERROR: no matching pattern found")
            elif len(spans) > 1:
                print(
                    "WARNING: there are multiple matches. Deleting before the first match")

        ind = spans[0][0]

        if verbose:
            if ind == 0:
                print(
                    "WARNING: the match already happens at the beginning -- not deleting anything.")

        if ind > 0:
            self.seq = self.seq[ind:]

    def extract_span(self, span):
        """Function to extract a substring given a span of indices."""
        substring = self.seq[span[0]:span[1] + 1]
        return substring

    def calculate_interface_affinity(self, span=None):
        """
        Function to calculate interface affinity of the sequence using the Wimley-White scale.
        Calculates delta_G of transfer from bulk water to POPC-water interface in kcal/mol.
        """
        if self.type != 'protein':
            print("ERROR: must be a protein sequence")
            raise TypeError
        if span is None:
            span = (0, len(self.seq))
        else:
            assert (len(span) == 2)
            span = (min(span), max(span))


        aff = self.get_value_array_based_on_rule(ww)
        int_affinity = [val * -1 for i, val in enumerate(aff) if i>=span[0] and i<=span[1]]

        return int_affinity

    def find_interface_affinity_period(self, window=5):
        """
        Function to find the periodicity in the protein interface affinity.
        """
        if self.type != 'protein':
            print("ERROR: must be a protein sequence")
            raise TypeError
        # since we calculate a moving average, any frequency above a certain freq will be meaningless
        int_affinity = self.calculate_interface_affinity()
        f, p = periodogram(movmean(int_affinity,window), 1)
        # exclude f=0
        f_res = f[np.argmax(p[1:])+1]
        return 1./f_res

    def calculate_weight(self):
        """Function to calculate the weight in Da of the given sequence."""
        assert (self.type == 'protein')
        weight = 0
        for aa in self.seq:
            weight += mw_aa[aa]

        return weight

    def calculate_composition(self):
        """Function to calculate the composition of a given sequence. Returns frequency of each letter (either amino acid or nucleic acid)."""
        if self.type == "protein":
            freqs = {i: 0 for i in aa_alphabet}
        elif self.type == "rna":
            freqs = {i: 0 for i in rna_alphabet}
        elif self.type == "dna":
            freqs = {i: 0 for i in dna_alphabet}
        else:
            print("ERROR: sequence type cannot be recognized")

        for aa in self.seq:
            freqs[aa] += 1

        return freqs

    def find_fragments(self, parent, frame, order, min_orf_len=90):
        """
        Given a sequence object (of type rna or dna), find orfs in a 
        specific reading frame.
        """
        seq_str = self.seq
        len_str = len(seq_str) // 3
        start_ind = 0
        end_ind = 0
        curr_fragment = ""
        recording = False
        orfs = []
        for ptr in range(0, 3*len_str, 3):
            curr_codon = seq_str[ptr:ptr+3]
            if curr_codon in start_codon and not recording:
                start_ind = ptr
                recording = True
            if recording:
                curr_fragment += curr_codon
            if curr_codon in stop_codon or ptr == 3*(len_str-1):
                if recording:
                    end_ind = ptr + 3
                    if len(curr_fragment) >= min_orf_len:
                        seq = sequence(curr_fragment, seq_type=self.type)
                        start_ind += frame
                        end_ind += frame
                        if order == -1:
                            start_ind = len(parent.seq) - start_ind
                            end_ind = len(parent.seq) - end_ind
                        if start_ind > end_ind:
                            start_ind, end_ind = end_ind, start_ind
                        orf = OpenReadingFrame(seq,
                                               parent,
                                               start_ind,
                                               end_ind,
                                               order,
                                               frame)
                        orfs.append(orf)
                    recording = False
                    curr_fragment = ""
        return orfs

    def detailed_orf_finder(self, min_orf_len=90):
        """
        Function to find open reading frames, and returns a list of OpenReadingFrame objects.
        """
        if self.type not in ('rna', 'dna'):
            print("ERROR: open reading frames can be found for DNA or RNA sequences only.")
            return None

        seq = sequence(self.seq.replace('-',''), seq_type=self.type)

        orfs = []
        for frame in range(3):
            for order in [1, -1]:
                if order == 1:
                    shifted = seq.frame_shift(frame=frame)
                else:
                    shifted = seq.reverse_complement().frame_shift(frame=frame)

                orfs_frame = shifted.find_fragments(seq, frame, order, min_orf_len=min_orf_len)
                for orf in orfs_frame:
                    orfs.append(orf)

        return orfs

    def find_index_after_alignment(self,
                                   unaligned_seq,
                                   index_without_gaps):
        """
        Function to find the index after alignment.
        """
        if not isinstance(unaligned_seq, sequence):
            raise TypeError
        if '-' in unaligned_seq.seq:
            print("ERROR: the unaligned sequence shouldn't have gaps.")
            return None
        
        seq_no_gaps = self.seq.replace('-','')
        start_idx = unaligned_seq.seq.find(seq_no_gaps)
        if start_idx == -1:
            print("ERROR: aligned and unaligned sequences do not match.")
            return None

        if index_without_gaps < start_idx:
            print("ERROR: the alignment doesn't contain the region of interest.")
            return None
        
        index_without_gaps -= start_idx
        
        return self.find_index_with_gaps(index_with_gaps=None,
                                            index_without_gaps=index_without_gaps)
        

    def find_index_with_gaps(self,
                             index_with_gaps=None,
                             index_without_gaps=None):
        """
        Given an index, convert to index without gaps or with gaps.
        """
        if index_with_gaps is None and index_without_gaps is None:
            print("ERROR: no index given.")
            return None

        if index_with_gaps is not None:
            # we need to convert an index with gaps to index without gaps
            if index_with_gaps >= len(self.seq):
                print(f"ERROR: index out of range for sequence of length {len(self.seq)}.")
                return None
            if self.seq[index_with_gaps] == '-':
                print("ERROR: given index corresponds to a gap.")
                return None
            substring = self.seq[:index_with_gaps+1]
            return len(substring.replace('-',''))-1

        if index_without_gaps is not None:
            # we need to convert an index without gaps to index with gaps
            seq_no_gaps = self.seq.replace('-','')
            if index_without_gaps >= len(seq_no_gaps):
                print(f"ERROR: index out of range for sequence without gaps of length {len(seq_no_gaps)}.")
                return None
            result = -1
            for i, char in enumerate(self.seq):
                if char != '-':
                    result += 1
                if result == index_without_gaps:
                    return i

    def find_start_codons(self):
        """
        Function that finds the start codons, assuming we are in the correct reading frame.
        """
        seq_str = self.seq
        indices = []
        n = len(seq_str)
        if self.type == 'protein':
            for ptr in range(n):
                if seq_str[ptr] == "M":
                    indices.append(ptr)
        elif self.type == 'dna' or self.type == 'rna':
            for ptr in range(0, n, 3):
                curr_codon = seq_str[ptr:ptr+3]
                if curr_codon in start_codon:
                    indices.append(ptr)
        else:
            print("ERROR: unknown sequence type.")
            return None

        return indices

    def find_stop_codons(self):
        """
        Function that finds the stop codons, assuming we are in the correct reading frame.
        """
        seq_str = self.seq
        indices = []
        n = len(seq_str)
        if self.type == 'protein':
            for ptr in range(n):
                if seq_str[ptr] == "M":
                    indices.append(ptr)
        elif self.type == 'dna' or self.type == 'rna':
            for ptr in range(0, n, 3):
                curr_codon = seq_str[ptr:ptr+3]
                if curr_codon in stop_codon:
                    indices.append(ptr)
        else:
            print("ERROR: unknown sequence type.")
            return None

        return indices

    def get_letter_frequencies(self):
        """
        Function that determines the frequencies of each letter (non-position-specific)
        and returns a dictionary.
        """
        if self.type == "protein":
            letters = {char:0 for char in aa_alphabet}
        elif self.type == "dna":
            letters = {char:0 for char in dna_alphabet}
        elif self.type == "rna":
            letters = {char:0 for char in rna_alphabet}
        else:
            print("WARNING: sequence type not recognized.")
            return None

        for char in self.seq:
            letters[char] += 1

        return letters

    def get_start_codon_nucleotide_frequencies(self, n_upstream=15, n_downstream=15, frames=None):
        """
        Function that determines start codons, and calculates frequencies of nucleotides
        at positions relative to the start codon.
        """
        if self.type != "rna" and self.type != "dna":
            print("ERROR: type has to be nucleotide")
            return None

        # init the frequencies dict
        frequencies = {}

        # make sure that the frames input is meaningful
        if frames is None:
            frames = range(3)
        elif isinstance(frames, list):
            # check that all we have is 0, 1, and 2
            if len(frames) > 3:
                print("ERROR: frames have to be a subset of [0, 1, 2].")
                raise ValueError
            for f in frames:
                if f not in [0, 1, 2]:
                    print("ERROR: the chosen frame has to be 0, 1, or 2.")
                    raise ValueError
        elif isinstance(frames, int):
            if frames not in [0, 1, 2]:
                print("ERROR: the chosen frame has to be 0, 1, or 2.")
                raise ValueError
            frames = [frames]
        else:
            print("ERROR: frames input not understood.")
            raise TypeError

        # process all requested frames
        for frame in frames:
            shifted = self.frame_shift(frame=frame)

            for ptr in range(0, len(shifted)+3, 3):
                curr_codon = shifted[ptr:ptr+3]
                if curr_codon in start_codon:
                    start_idx = ptr - n_upstream
                    end_idx = ptr + n_downstream + 3
                    if start_idx >= 0 and end_idx < len(shifted):
                        selected = shifted[start_idx:end_idx]
                        for i, char in enumerate(selected):
                            pos = i - n_upstream
                            key = (char, pos)
                            if key in frequencies.keys():
                                frequencies[key] += 1
                            else:
                                frequencies[key] = 1
        return frequencies

    def get_start_codon_nucleotide_frequencies_re(self, n_upstream=15, n_downstream=15):
        """
        The same as above, but with regex.
        """

        # compile re
        if self.type == 'rna':
            p = re.compile('AUG')
            possible_chars = ['A','U','C','G','N']
        else:
            p = re.compile('ATG')
            possible_chars = ['A','T','C','G','N']
        seq_str = self.seq

        # initialize dict
        frequencies = {}
        positions = list(range(-1*n_upstream,0,1)) + list(range(0,n_downstream+4))
        for char in possible_chars:
            for pos in positions:
                frequencies[(char, pos)] = 0

        for match in p.finditer(seq_str):
            start_idx = match.start() - n_upstream
            end_idx = match.end() + n_downstream 
            if start_idx >= 0 and end_idx < len(seq_str):
                selected = seq_str[start_idx:end_idx]
                for i, char in enumerate(selected):
                    pos = i - n_upstream
                    key = (char, pos)
                    frequencies[key] += 1

        return frequencies

    def write_fasta(self, file_name, append=True):
        """
        Function that writes sequence to a .fasta file. If append
        is True, the .fasta file is appended by the sequence.

        Bypasses the creation of seq_set to write sequences to file.
        Useful for pyspark and other large set manipulations.
        """
        if append:
            file = open(file_name, 'a', encoding='utf-8')
        else:
            file = open(file_name, 'w', encoding='utf-8')

        if self.name is not None:
            name = ">" + self.name
        else:
            name = ">"

        file.write(name)
        file.write('\n')
        write_in_columns(file, self.seq, ncols=79)
        file.close()

    def find_indices_of_match(self, regex):
        """
        Given a regular expression, find all indices of match and return.
        (index of the beginning of the match)
        """
        if isinstance(regex, str):
            pattern = re.compile(regex)
        elif isinstance(regex, re.Pattern):
            pattern = regex
        else:
            print(f"ERROR: input regex of {type(regex)} not supported.")
            raise TypeError
        indices = []
        for x in pattern.finditer(self.seq):
            indices.append(x.span()[0])

        return indices

    def get_value_array_based_on_rule(self, rule_dict):
        """
        Given a dictionary of rules that assigns values to letters,
        return an array of values that correspond to the sequence.
        """
        keys = set(list(rule_dict.keys()))
        wrong = False
        if self.type == 'protein':
            if keys != set(aa_alphabet):
                wrong = True
        elif self.type == 'dna':
            if keys != set(dna_alphabet):
                wrong = True
        elif self.type == 'rna':
            if keys != set(rna_alphabet):
                wrong = True
        else:
            print(f"ERROR: wrong sequence type ({self.type}).")
            raise TypeError
        if wrong:
            print(f"ERROR: rule dictionary doesn't have keys consistent with sequence type {self.type}.")
            raise KeyError

        vals = []
        for char in self.seq:
            vals.append(rule_dict[char])

        return vals


    def find_period_based_on_rule(self, rule_dict, period_range=(5,200), window=5):
        """
        Given a rule dictionary that assigns values to letters, return
        the most dominant period from periodogram.

        Peaks in periodogram within the given period range are considered.
        """
        if not isinstance(rule_dict, dict):
            print("ERROR: rule_dict should be of type dict.")
            raise TypeError
        vec = self.get_value_array_based_on_rule(rule_dict=rule_dict)

        freq_range = [1./max(period_range), 1./min(period_range)]

        f, p = periodogram(movmean(vec, window), 1)
        p_max = 0
        f_max = 0
        for i, f_val in enumerate(f):
            if f_val >= freq_range[0] and f_val <= freq_range[1]:
                if p[i] > p_max:
                    p_max = p[i]
                    f_max = f_val

        return 1./f_max

    def find_kmers(self, k):
        """
        Function to find kmers in the sequence with their frequencies.
        """
        return find_kmers_in_string(self.seq, k)

def read_single_seq_from_file(file_name, seq_name=None, seq_idx=None):
    """
    Given a large .fasta file with many sequences, read and return a single
    sequence only.
    """
    if seq_name is None and seq_idx is None:
        print("ERROR: you need to provide either a sequence name or an index.")
        raise ValueError

    if seq_name is not None and seq_idx is not None:
        print("WARNING: both a sequence name and index has been provided. I will use the name and ignore the index.")
        seq_idx = None

    count = 0
    seq_str = ""
    started = False
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith(">"):
                # we have the beginning of a sequence
                if started:
                    break
                if (seq_name is not None and seq_name in line) or count == seq_idx:
                    name = line.strip()[1:]
                    started = True
                count += 1
            else:
                if started:
                    seq_str += line.strip()

    if not started or seq_str == "":
        print("WARNING: no match found.")
        return None

    seq = sequence(seq_str, name=name)
    return seq