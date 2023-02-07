"""
This file holds the sequence class and related methods.
"""
from Bio.Seq import Seq
import sequence_analysis.utils as utils
import re
from sequence_analysis.utils import aa_alphabet
from sequence_analysis.utils import dna_alphabet
from sequence_analysis.utils import rna_alphabet
from sequence_analysis.utils import diff_letters
from colorama import Fore, Back, Style
from sequence_analysis.utils import ww
from  sequence_analysis.utils import mw_aa

class sequence:
    # TODO: inherit from the Seq class?
    def __init__(self,seq,name=None,type=None):
        self.seq = seq.upper()
        self.name = name
        self.type = type

        if type is None:
            self.set_type()

    def __eq__(self, other):
        if type(self) == type(other):
            if self.type == other.type:
                return self.seq == other.seq
            else:
                return False
        else:
            return False

    def __gt__(self, other):
        assert(type(self) == type(other))
        assert(self.type == other.type)

        if self.seq > other.seq:
            return True
        else:
            return False

    def __lt__(self, other):
        assert(type(self) == type(other))
        assert(self.type == other.type)

        if self.seq < other.seq:
            return True
        else:
            return False

    def __str__(self):
        return f"{self.seq}"

    def __repr__(self):
        return f"<Sequence object of type {self.type} with '{self.seq:.5}...' at {hex(id(self))}>" 

    def set_type(self):
        s = self.seq.upper()
        all_letters = set([*s])

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


    def frame_shift(self, frame=0):
        if self.type != 'rna' and self.type != 'dna':
            print(f"ERROR: type {self.type} is not meaningful for frame shifting. Please provide RNA or DNA sequences.")
            return None

        seq = self.seq

        # returns frame shifted sequence
        if frame == 0:
            if len(seq) % 3 == 0:
                return seq
        elif frame == 1:
            seq = seq[1:]
        elif frame == 2:
            seq = seq[2:]
        else:
            print("ERROR: problem with frame shifting. Frame should be 0, 1, or 2.")

        N_to_be_added = 3 - len(seq) % 3
        seq = seq + "N" * N_to_be_added

        # TODO: move to the test?
        assert(len(seq) % 3 == 0)

        return seq

    def six_frame_check(self, regex):
        # does a 6-frame translation and checks for the existence of a given regex pattern
        selected = False
        for frame in range(3):
            for order in range(1):
                if order == 0:
                    s = Seq(self.frame_shift(frame=frame)).translate()
                    print(s)
                else:
                    s = Seq(self.frame_shift(frame=frame)).reverse_complement().translate()
                    print(s)

                # create sequence object with the translated sequence
                seq = sequence(str(s))
                print(frame,order,seq.seq)

                # NOTE: the following implicitly assumes that only one reading frame is valid
                selected = seq.check_for_pattern(regex)
                if selected:
                    true_sequence = str(s)
                    return true_sequence

        return None

    def check_for_pattern(self, regex):
        # TODO: check if regex and seq from the same alphabet?
        return utils.check_for_pattern(self.seq, regex)

    def choose_all_matching_patterns(self, regex, return_between_matching=False):
        s = self.seq
        p = re.compile(regex)

        x = p.findall(s)

        spans = utils.find_spans(s,x)

        if return_between_matching:
            inverted = utils.invert_spans(spans)
            return x, spans, inverted        

        return x, spans

    def choose_all_matching_patterns_and_print(self, regex):
        matching, spans, inverted = self.choose_all_matching_patterns(regex,return_between_matching=True)

        seq = self.seq

        to_print = f""
        initial = seq[:spans[0][0]]
        to_print += initial
        for i, span in enumerate(spans):
            to_print += f"{Fore.BLUE}" + seq[span[0]:span[1]+1] + f"{Style.RESET_ALL}"
            if i < len(spans) - 1:
                to_print += f"{Fore.RED}" + seq[inverted[i][0]:inverted[i][1]+1] + f"{Style.RESET_ALL}"
        tail = seq[spans[-1][1]+1:]
        to_print += tail

        print(to_print)

    def remove_before_pattern(self, regex, verbose=True):
        x, spans = self.choose_all_matching_patterns(regex)
        if verbose:
            if len(spans) == 0:
                print("ERROR: no matching pattern found")
            elif len(spans) > 1:
                print("WARNING: there are multiple matches. Deleting before the first match")
        
        ind = spans[0][0] 

        if verbose:
            if ind == 0:
                print("WARNING: the match already happens at the beginning -- not deleting anything.")

        if ind > 0:
            self.seq = self.seq[ind:]


    def extract_span(self,span):
        substring = self.seq[span[0]:span[1]+1]
        return substring

    def calculate_interface_affinity(self, span=None, kcal_per_mol=True):
        if span is None:
            span = (0, len(self.seq))
        else:
            assert(len(span)==2)
            span = (min(span), max(span))

        assert(self.type == 'protein')

        seq = []
        int_affinity = []

        for i in range(span[0],span[1]):
            aa = self.seq[i]
            int_affinity.append(ww[aa]*-1)
            seq.append(aa)

        return int_affinity, seq

    def calculate_weight(self):
        assert(self.type == 'protein')
        weight = 0
        for aa in self.seq:
            weight += mw_aa[aa]

        return weight

    def calculate_composition(self):
        if self.type == "protein":
            freqs = {i:0 for i in aa_alphabet}
        elif self.type == "rna":
            freqs = {i:0 for i in rna_alphabet}
        elif self.type == "dna":
            freqs = {i:0 for i in dna_alphabet}
        else:
            print("ERROR: sequence type cannot be recognized")

        for aa in self.seq:
            freqs[aa] += 1

        return freqs