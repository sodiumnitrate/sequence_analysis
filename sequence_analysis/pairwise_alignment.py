from sequence_analysis.utils import query_blosum50

class pairwise_alignment:
    def __init__(self, sequence1, sequence2, algorithm="needleman-wunsch"):
        # TODO: check input data structures
        self.sequence1 = sequence1
        self.sequence2 = sequence2

        self.algorithm = algorithm

        self.result_string = None
        self.score = None
        self.F = None

    def align(self):
        if self.algorithm == "needleman-wunsch":
            self.needleman_wunsch()

    def needleman_wunsch(self,d=-8):
        # get lengths of sequences
        n = len(self.sequence1.seq)
        m = len(self.sequence2.seq)

        # initialize matrix
        F = [[0 for _ in range(n)] for _ in range(m)]
        pointers = [["" for _ in range(n)] for _ in range(m)]

        # initialize the boundaries
        for i in range(n):
            F[i][0] = -i * d

        for i in range(m):
            F[0][i] = -i * d
