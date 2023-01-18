from sequence_analysis.utils import query_blosum50
from sequence_analysis.sequence import sequence
import pdb

class pairwise_alignment:
    def __init__(self, sequence1, sequence2, algorithm="needleman-wunsch"):
        # TODO: check input data structures
        if isinstance(sequence1, str):
            sequence1 = sequence(sequence1)
        if isinstance(sequence2, str):
            sequence2 = sequence(sequence2)
        self.sequence1 = sequence1
        self.sequence2 = sequence2

        self.algorithm = algorithm

        self.sequence1_aligned = None
        self.sequence2_aligned = None
        self.score = None
        self.F = None
        self.pointers = None

    def align(self,d=8,verbose=False):
        if self.algorithm == "needleman-wunsch":
            self.needleman_wunsch(d=d,verbose=verbose)
        elif self.algorithm == "smith_waterman":
            self.smith_waterman(d=d, verbose=verbose)

    def needleman_wunsch(self,d=8,verbose=False):
        # get lengths of sequences
        n = len(self.sequence1.seq)
        m = len(self.sequence2.seq)

        # initialize matrix
        F = [[0 for _ in range(n+1)] for _ in range(m+1)]
        pointers = [["" for _ in range(n+1)] for _ in range(m+1)]

        # initialize the boundaries
        for i in range(m+1):
            F[i][0] = -i * d
            pointers[i][0] = 'u'

        for i in range(n+1):
            F[0][i] = -i * d
            pointers[0][i] = 'l'

        pointers[0][0] = ""

        # TODO: rewrite this in C++ because nested for loops
        # i indexes into y, j indexes into x
        for i in range(1,m+1):
            for j in range(1,n+1):
                possible_values = []

                # x_j and y_i are aligned
                aa1 = self.sequence1.seq[j-1]
                aa2 = self.sequence2.seq[i-1]
                possible_values.append(F[i-1][j-1] + query_blosum50(aa1,aa2))

                # x_i is aligned to a gap
                possible_values.append(F[i-1][j] - d)

                # y_j is aligned to a gap
                possible_values.append(F[i][j-1] - d)

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
            self.print_matrix(type="pointers")

        self.traceback_nw()

    def traceback_nw(self):
        n = len(self.sequence1.seq)
        m = len(self.sequence2.seq)

        x_str = ""
        y_str = ""

        i = m
        j = n
        while i > 0 or j > 0:
            # the boundary conditions on F should prevent from either index from going to negative values
            direction = self.pointers[i][j]
            if direction == "u":
                # y_j is aligned with a gap
                x_str = "-" + x_str
                y_str = self.sequence2.seq[i-1] + y_str
                i -= 1
            elif direction == "l":
                # x_i is aligned with a gap
                x_str = self.sequence1.seq[j-1] + x_str
                y_str = "-" + y_str
                j -= 1
            elif direction == "ul":
                # x_i is aligned with y_j
                x_str = self.sequence1.seq[j-1] + x_str
                y_str = self.sequence2.seq[i-1] + y_str
                i -= 1
                j -= 1

        self.sequence1_aligned = sequence(x_str)
        self.sequence2_aligned = sequence(y_str)

    def print_matrix(self,type="score"):
        # method to print either the score or pointer matrix for debugging/visualization purposes
        if type == "score":
            matrix = self.F
        elif type == "pointers":
            matrix = self.pointers
        else:
            print("ERROR: please specify 'score' or 'pointers'")
            return 

        if matrix is None:
            print("ERROR: matrix is empty. Did you run align() first?")
            return

        for line in matrix:
            print(line)

    def smith_waterman(self,d=8,verbose=False):
        n = len(self.sequence1.seq)
        m = len(self.sequence2.seq)

        # initialize matrix
        F = [[0 for _ in range(m+1)] for _ in range(n+1)]
        pointers = [["" for _ in range(m+1)] for _ in range(n+1)]

        # TODO: rewrite in C++
        for i in range(1,n+1):
            for j in range(1,m+1):
                possible_values = []
                aa1 = self.sequence1.seq[i-1]
                aa2 = self.sequence2.seq[j-1]

                # start a new alignment
                possible_values.append(0)

                # x_i and y_j are aligned
                possible_values.append(F[i-1][j-1] + query_blosum50(aa1,aa2))

                # x_i is aligned to gap
                possible_values.append(F[i-1][j] - d)

                # y_j is aligned to gap
                possible_values.append(F[i][j-1] - d)

                # choose option that maximizes the cumulative score
                max_val = max(possible_values)
                max_ind = possible_values.index(max_val)

                F[i][j] = max_val
                if max_ind == 0:
                    pointers[i][j] == ""
                elif max_ind == 1:
                    pointers[i][j] == "ul"
                elif max_ind == 2:
                    pointers[i][j] == "l"
                elif max_ind == 3:
                    pointers[i][j] == "u"

        self.score = max(max(F, key=max))
        self.F = F
        self.pointers = pointers

        if verbose:
            self.print_matrix()
            self.print_matrix(type="pointers")
        self.traceback_sw()

    def traceback_sw(self):
        # find index of max score
        pass