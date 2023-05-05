"""
This module holds the pairwise_alignment object and its relevant methods,
most importantly the align() method. The most efficient global alignment
setting uses biopython's pairwise alignment methods under the hood.

Note that my own (purely python) implementations of Needleman-Wunsch and
Smith-Waterman are also in this module for educational purposes.
"""
from unicodedata import name
from Bio import Align
from sequence_analysis.utils import query_blosum50
from sequence_analysis.sequence import sequence
from sequence_analysis.utils import blosum_50, blosum_62


class pairwise_alignment:
    """
    This class sets up pairwise alignment with the given parameters.
    After the necessary parameters have been set up, the align() method
    can be run to compute pairwise alignments.
    """

    def __init__(self, target, query, match=1, unmatch=0, gap=-8,
                 gap_open=-9, use_blosum_50=False, algorithm="needleman-wunsch",
                 use_blosum_62=False):
        if isinstance(target, str):
            target = sequence(target)
        elif not isinstance(target, sequence):
            raise TypeError
        if isinstance(query, str):
            query = sequence(query)
        elif not isinstance(query, sequence):
            raise TypeError
        
        self.target = target
        self.query = query

        # what parameters to use.
        # if use_blosum_50 is True, the match and unmatch params are ignored
        self.use_blosum_50 = use_blosum_50
        self.use_blosum_62 = use_blosum_62
        self.match = match
        self.unmatch = unmatch
        # gap is -8 by default, which is chosen to play nicely with blosum50.
        # if using custom match and unmatch scores, make sure the chosen gap
        # penalty is reasonable
        self.gap = gap
        self.gap_open = gap_open

        self.algorithm = algorithm

        self.target_aligned = None
        self.query_aligned = None
        self.score = None
        self.F = None
        self.pointers = None

    def check_input_params(self):
        """
        This function makes sure that the user-input parameters
        are reasonable.
        """
        # make sure both sequences have the same type
        type_target = self.target.type
        type_query = self.query.type

        if type_target is None or type_query is None:
            raise TypeError

        if type_target != type_query:
            raise TypeError

        # make sure blosum is requested for protein only
        protein = type_target == 'protein'
        if not protein and (self.use_blosum_50 or self.use_blosum_62):
            raise TypeError

        # make sure both blosum matrices aren't requested
        if self.use_blosum_50 and self.use_blosum_62:
            print("WARNING: BLOSUM50 and BLOSUM62 are both requested. Using the latter.")
            self.use_blosum_50 = False

    def align(self, verbose=False):
        """Function that aligns the two sequences held in pairwise_alignment."""

        # first, check that the input params make sense
        self.check_input_params()

        if self.algorithm == "needleman-wunsch":
            self.needleman_wunsch(verbose=verbose)
        elif self.algorithm == "smith-waterman":
            self.smith_waterman(verbose=verbose)
        elif self.algorithm == "biopython-global":
            self.biopython_global()
        elif self.algorithm == "biopython-local":
            self.biopython_local()
        else:
            print(
                f"ERROR: algorithm {self.algorithm} is not recognized. Not aligning.")

    def biopython_global(self):
        """Function that sets up a global alignment usin biopython's PairwiseAligner."""
        target = self.target.seq
        query = self.query.seq

        # setup aligner
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'

        if self.use_blosum_50:
            aligner.substitution_matrix = blosum_50
        elif self.use_blosum_50:
            aligner.substitution_matrix = blosum_62
        else:
            aligner.match_score = self.match
            aligner.mismatch_score = self.unmatch
            aligner.open_gap_score = self.gap_open
            aligner.extend_gap_score = self.gap

        alignments = aligner.align(target, query)
        self.score = alignments.score
        self.target_aligned = sequence(alignments[0][0],
                                       seq_type=self.target.type,
                                       name=self.target.name)
        self.query_aligned = sequence(alignments[0][1],
                                      seq_type=self.query.type,
                                      name=self.target.name)

    def biopython_local(self):
        """Function that sets up a local alignment using biopython's PairwiseAligner."""
        seq1 = self.target.seq
        seq2 = self.query.seq

        # setup
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'

        if self.use_blosum_50:
            aligner.substitution_matrix = blosum_50
        elif self.use_blosum_62:
            aligner.substitution_matrix = blosum_62
        else:
            aligner.match_score = self.match
            aligner.mismatch_score = self.unmatch
            aligner.open_gap_score = self.gap_open
            aligner.extend_gap_score = self.gap

        alignments = aligner.align(seq1, seq2)

        self.score = alignments.score
        self.target_aligned = sequence(alignments[0][0],
                                       seq_type=self.target.type,
                                       name=self.target.name)
        self.query_aligned = sequence(alignments[0][1],
                                      seq_type=self.query.type,
                                      name=self.target.name)

    def needleman_wunsch(self, verbose=False):
        """
        Function that calculates the score matrix for Needleman-Wunsch.

        Warning: exists for educational purposes. Significantly slower
        than biopython-global.
        """
        # get lengths of sequences
        n = len(self.target.seq)
        m = len(self.query.seq)

        # initialize matrix
        F = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
        pointers = [["" for _ in range(n + 1)] for _ in range(m + 1)]

        # initialize the boundaries
        for i in range(m + 1):
            F[i][0] = i * self.gap
            pointers[i][0] = 'u'

        for i in range(n + 1):
            F[0][i] = i * self.gap
            pointers[0][i] = 'l'

        pointers[0][0] = ""

        # TODO: rewrite this in C++ because nested for loops
        # i indexes into y, j indexes into x
        if self.use_blosum_50:
            # TODO: REFACTOR
            for i in range(1, m + 1):
                for j in range(1, n + 1):
                    possible_values = []

                    # x_j and y_i are aligned
                    aa1 = self.target.seq[j - 1]
                    aa2 = self.query.seq[i - 1]
                    possible_values.append(
                        F[i - 1][j - 1] + query_blosum50(aa1, aa2))

                    # x_i is aligned to a gap
                    possible_values.append(F[i - 1][j] + self.gap)

                    # y_j is aligned to a gap
                    possible_values.append(F[i][j - 1] + self.gap)

                    # choose option that maximizes the cumulative score
                    max_val = max(possible_values)
                    max_ind = possible_values.index(max_val)

                    F[i][j] = max_val
                    if max_ind == 0:
                        pointers[i][j] = "ul"
                    elif max_ind == 1:
                        pointers[i][j] = "u"
                    else:
                        pointers[i][j] = "l"
        else:
            for i in range(1, m + 1):
                for j in range(1, n + 1):
                    possible_values = []

                    # x_j and y_i are aligned
                    aa1 = self.target.seq[j - 1]
                    aa2 = self.query.seq[i - 1]
                    if aa1 == aa2:
                        sc = self.match
                    else:
                        sc = self.unmatch
                    possible_values.append(F[i - 1][j - 1] + sc)

                    # x_i is aligned to a gap
                    possible_values.append(F[i - 1][j] + self.gap)

                    # y_j is aligned to a gap
                    possible_values.append(F[i][j - 1] + self.gap)

                    # choose option that maximizes the cumulative score
                    max_val = max(possible_values)
                    max_ind = possible_values.index(max_val)

                    F[i][j] = max_val
                    if max_ind == 0:
                        pointers[i][j] = "ul"
                    elif max_ind == 1:
                        pointers[i][j] = "u"
                    else:
                        pointers[i][j] = "l"

        self.score = F[-1][-1]
        self.F = F
        self.pointers = pointers

        if verbose:
            self.print_matrix()
            self.print_matrix(info_type="pointers")

        self.traceback_nw()

    def traceback_nw(self):
        """
        Function that performs traceback on the score matrix of Needleman-Wunsch to
        obtain the final alignment strings.
        """
        n = len(self.target.seq)
        m = len(self.query.seq)

        x_str = ""
        y_str = ""

        i = m
        j = n
        while i > 0 or j > 0:
            # the boundary conditions on F should prevent from either index
            # from going to negative values
            direction = self.pointers[i][j]
            if direction == "u":
                # y_j is aligned with a gap
                x_str = "-" + x_str
                y_str = self.query.seq[i - 1] + y_str
                i -= 1
            elif direction == "l":
                # x_i is aligned with a gap
                x_str = self.target.seq[j - 1] + x_str
                y_str = "-" + y_str
                j -= 1
            elif direction == "ul":
                # x_i is aligned with y_j
                x_str = self.target.seq[j - 1] + x_str
                y_str = self.query.seq[i - 1] + y_str
                i -= 1
                j -= 1

        self.target_aligned = sequence(x_str,
                                       seq_type=self.target.type,
                                       name=self.target.name)
        self.query_aligned = sequence(y_str,
                                      seq_type=self.query.type,
                                      name=self.query.name)

    def print_matrix(self, info_type="score"):
        """Function that prints either the score or pointer matrix."""
        # method to print either the score or pointer matrix for
        # debugging/visualization purposes
        if info_type == "score":
            matrix = self.F
        elif info_type == "pointers":
            matrix = self.pointers
        else:
            print("ERROR: please specify 'score' or 'pointers'")
            return

        if matrix is None:
            print("ERROR: matrix is empty. Did you run align() first? Note that the biopython implementation does not return this matrix.")
            return

        for line in matrix:
            print(line)

    def smith_waterman(self, verbose=False):
        """Function that calculates the score matrix for local alignment with Smith-Waterman."""
        # TODO: get rid of this properly
        print("WARNING: this function is not tested for custom alignment scoring.")
        n = len(self.target.seq)
        m = len(self.query.seq)

        # initialize matrix
        F = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
        pointers = [["" for _ in range(n + 1)] for _ in range(m + 1)]

        # TODO: rewrite in C++
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                possible_values = []
                aa1 = self.target.seq[j - 1]
                aa2 = self.query.seq[i - 1]

                # start a new alignment
                possible_values.append(0)

                # x_i and y_j are aligned
                possible_values.append(
                    F[i - 1][j - 1] + query_blosum50(aa1, aa2))

                # x_i is aligned to gap
                possible_values.append(F[i - 1][j] + self.gap)

                # y_j is aligned to gap
                possible_values.append(F[i][j - 1] + self.gap)

                # choose option that maximizes the cumulative score
                max_val = max(possible_values)
                max_ind = possible_values.index(max_val)

                F[i][j] = max_val
                if max_ind == 0:
                    pointers[i][j] = ""
                elif max_ind == 1:
                    pointers[i][j] = "ul"
                elif max_ind == 2:
                    pointers[i][j] = "u"
                elif max_ind == 3:
                    pointers[i][j] = "l"

        self.score = max(max(F, key=max))
        self.F = F
        self.pointers = pointers

        if verbose:
            self.print_matrix()
            self.print_matrix(info_type="pointers")
        self.traceback_sw()

    def traceback_sw(self):
        """
        Function that performs traceback on the score matrix of Smith-Waterman to
        obtain the final alignment strings.
        """
        i_start = None
        j_start = None
        # find index of max score
        for i in range(len(self.F)):
            for j in range(len(self.F[0])):
                if self.F[i][j] == self.score:
                    i_start = i
                    j_start = j
                    break
            if i_start is not None:
                break

        assert (i_start is not None)

        x_str = ""
        y_str = ""

        curr_score = self.F[i_start][j_start]
        while curr_score != 0:
            direction = self.pointers[i_start][j_start]
            if direction == "u":
                x_str = "-" + x_str
                y_str = self.query.seq[i_start - 1] + y_str
                i_start -= 1
            elif direction == "l":
                x_str = self.target.seq[j_start - 1] + x_str
                y_str = "-" + y_str
                j_start -= 1
            elif direction == "ul":
                x_str = self.target.seq[j_start - 1] + x_str
                y_str = self.query.seq[i_start - 1] + y_str
                i_start -= 1
                j_start -= 1

            curr_score = self.F[i_start][j_start]

        self.target_aligned = sequence(x_str,
                                       seq_type=self.target.type,
                                       name=self.target.name)
        self.query_aligned = sequence(y_str,
                                      seq_type=self.query.type,
                                      name=self.query.name)
