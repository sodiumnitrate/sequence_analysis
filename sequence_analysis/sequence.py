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

    def choose_all_matching_patterns(self, regex):
        s = self.seq
        p = re.compile(regex)

        x = p.findall(s)

        '''
        spans = []
        for match in x:
            # TODO: edit this to make sure duplicates are treated correctly
            # regex to match the match
            p2 = re.compile(match)
            # search seq for match
            x2 = p2.search(s)
            # get range of chards for match
            span = x2.span()
            # append to list of spans
            spans.append(span)

        # sort spans based on their lower bound
        spans.sort(key=lambda y: y[1])

        # merge any overlapping spans (not sure why this is happening in the first place)
        new_spans = []
        for span1 in spans:
            overlap = False
            for span2 in spans:
                # check if the two spans overlap
                cond1 = (span1[0] < span2[0] and span1[1] > span2[0]) 
                cond2 = (span2[0] < span1[0] and span2[1] > span1[0])
                if cond1 or cond2:
                    merged = (min(span1[0],span2[0]), max(span1[1],span2[1]))
                    overlap = True
                    if merged not in new_spans:
                        new_spans.append(merged)
            if not overlap:
                new_spans.append(span1)
        '''

        spans = utils.find_spans(s,x)

        return x, spans
