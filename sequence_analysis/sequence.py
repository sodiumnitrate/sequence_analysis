"""
This file holds the sequence class and related methods.
"""
from Bio.Seq import Seq
import sequence_analysis.utils as utils
import re

class sequence:
    # TODO: inherit from the Seq class?
    def __init__(self,seq,name=None,type=None):
        self.seq = seq.upper()
        self.name = name
        self.type = type

    def frame_shift(self, frame=0):
        if self.type != 'rna' and self.type != 'dna':
            print(f"ERROR: type {self.type} is not meaningful for frame shifting. Please provide RNA or DNA sequences.")
            return None

        seq = self.seq

        # returns frame shifted sequence
        if frame == 0:
            pass

        if frame == 1:
            seq = seq[1:]

        if frame == 2:
            seq = seq[2:]

        N_to_be_added = 3 - len(seq) % 3
        seq = seq + "N" * N_to_be_added

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

                # NOTE: the following implicitly assumes that only one reading frame is valid
                selected = seq.check_for_pattern(regex)
                if selected:
                    true_sequence = str(s)
                    return true_sequence

        return None

    def check_for_pattern(self, regex):
        # TODO: rethink this wrapper situation
        # TODO: check if regex and seq from the same alphabet?
        return utils.check_for_pattern(self.seq, regex)
