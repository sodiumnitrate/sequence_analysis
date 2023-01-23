"""
This file holds the sequence class and related methods.
"""
from Bio.Seq import Seq
import sequence_analysis.utils as utils
import re
from sequence_analysis.utils import dna_alphabet
from sequence_analysis.utils import rna_alphabet
from sequence_analysis.utils import diff_letters

class sequence:
    # TODO: inherit from the Seq class?
    def __init__(self,seq,name=None,type=None):
        self.seq = seq.upper()
        self.name = name
        self.type = type

        if type is None:
            self.set_type()

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


