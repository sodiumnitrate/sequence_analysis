"""
This file holds the sequence class and related methods.
"""

from typing import Sequence
from Bio.Seq import Seq
import sequence_analysis.utils as utils
import re
from sequence_analysis.utils import aa_alphabet
from sequence_analysis.utils import dna_alphabet
from sequence_analysis.utils import rna_alphabet
from sequence_analysis.utils import diff_letters
from colorama import Fore, Back, Style
from sequence_analysis.utils import ww
from sequence_analysis.utils import mw_aa
from sequence_analysis.utils import start_codon
from sequence_analysis.utils import stop_codon
from sequence_analysis.utils import stop_codon_re
from sequence_analysis.open_reading_frame import OpenReadingFrame

class sequence:
    """This class corresponds to a single biological sequence (protein, rna, or dna)"""
    # TODO: implement slicing

    def __init__(self, seq, name=None, seq_type=None):
        """This function initializes the sequence object."""
        self.seq = seq.upper()
        self.name = name
        self.type = seq_type

        if seq_type is None:
            self.set_type()

    def __eq__(self, other):
        """
        This function overloads the __eq__ function such that if two sequences are of the
        same type and have the same exact string for sequence, they are deemed equal.
        """
        if isinstance(self, type(other)):
            if self.type == other.type:
                return self.seq == other.seq
            else:
                return False
        else:
            return False

    def __gt__(self, other):
        """This function overloads __gt__ to use the string __gt__ on the seq attribute."""
        assert (isinstance(self, type(other)))

        if self.type != other.type:
            if not(self.type is None and other.type is None):
                return False

        if self.seq > other.seq:
            return True
        else:
            return False

    def __lt__(self, other):
        """"This function overloads __lt__ to use the string __lt__ on the seq attribute."""
        assert (isinstance(self, type(other)))

        if self.type != other.type:
            if not(self.type is None and other.type is None):
                return False

        if self.seq < other.seq:
            return True
        else:
            return False

    def __str__(self):
        """Function that overloads __str__ for sequence object."""
        return f"{self.seq}"

    def __repr__(self):
        """Function that overloas __repr__ for sequence object."""
        return f"<Sequence object of type {self.type} with '{self.seq:.5}...' at {hex(id(self))}>"

    def __len__(self):
        """Overwrites __len__ to return the number of chars in seq.seq"""
        return len(self.seq)

    def set_type(self, seq_type=None):
        """Function that sets the type of sequence based on the letters it contains."""
        # TODO: add custom input option
        # TODO: add "N" to the rna and dna alphabet?
        s = self.seq.upper()
        if len(s) == 0:
            print("WARNING: empty sequence. Can't set type.")
            self.type = None
            return
        all_letters = set([*s])

        if seq_type is not None:
            self.type = seq_type
            return

        if len(diff_letters.intersection(all_letters)):
            self.type = 'protein'
        elif set(dna_alphabet).issubset(all_letters) or set(rna_alphabet).issubset(all_letters):
            if 'U' in all_letters:
                self.type = 'rna'
            else:
                # this could fail if we have RNA sequences that happens to have
                # no Us, but that's unlikely
                self.type = 'dna'
        else:
            print("I can't assign a type, please specify manually.")

    def translate(self):
        """
        Function that translates DNA or RNA sequence into protein.
        """
        if self.type != 'rna' and self.type != 'dna':
            print(f"WARNING: no translation for sequence type {self.type}.")
            return
        s = Seq(self.seq.replace('-',""))
        s = str(s.translate())
        s = s.strip("X")

        return sequence(s)

    def find_codon(self, codon):
        """
        Function to search for codons. Returns index of first match.
        Returns None if not rna or dna.
        Returns -1 if no match.
        """
        if self.type != "rna" and self.type != "dna":
            print("ERROR: cannot search for codons if the sequence isn't rna or dna.")
            return None

        seq_str = self.seq
        codons = [seq_str[3*ptr:3*ptr+3] for ptr in range(len(seq_str) // 3)]
        try:
            ind = codons.index(codon)
        except ValueError:
            return -1

        return ind * 3

    def reverse_complement(self):
        if self.type != "rna" and self.type != "dna":
            print(f"ERROR: cannot reverse complement: the sequence is of type {self.type}.")
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
        return sequence(seq, seq_type=self.type)

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

    def check_for_pattern(self, regex):
        """Function that checks for a given regex pattern in the sequence."""
        # TODO: check if regex and seq from the same alphabet?
        return utils.check_for_pattern(self.seq, regex)

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
                to_print += f"{Fore.RED}" + seq[inverted[i][0]                                                :inverted[i][1] + 1] + f"{Style.RESET_ALL}"
        tail = seq[spans[-1][1] + 1:]
        to_print += tail

        print(to_print)

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

    def calculate_interface_affinity(self, span=None, kcal_per_mol=True):
        """Function to calculate interface affinity of the sequence using the Wimley-White scale."""
        if span is None:
            span = (0, len(self.seq))
        else:
            assert (len(span) == 2)
            span = (min(span), max(span))

        assert (self.type == 'protein')

        seq = []
        int_affinity = []

        for i in range(span[0], span[1]):
            aa = self.seq[i]
            int_affinity.append(ww[aa] * -1)
            seq.append(aa)

        return int_affinity, seq

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
            if curr_codon in stop_codon:
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

        orfs = []
        for frame in range(3):
            for order in [1, -1]:
                if order == 1:
                    shifted = self.frame_shift(frame=frame)
                else:
                    shifted = self.reverse_complement().frame_shift(frame=frame)

                orfs_frame = shifted.find_fragments(self, frame, order, min_orf_len=min_orf_len)
                for orf in orfs_frame:
                    orfs.append(orf)

        return orfs