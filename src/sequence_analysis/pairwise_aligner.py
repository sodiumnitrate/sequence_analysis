"""
This file holds the PairwiseAligner class and related methods.
"""
from .sequence_analysis_cpp import PairwiseAligner as PairwiseAligner_cpp

class PairwiseAligner(PairwiseAligner_cpp):
    """
    Python bindings for the PairwiseAligner class.

    std::string algorithm
    std::string target
    std::string query
    float score
    """
    def __init__(self, scoring="blosum50", gap_penalty=-8, target=None, query=None, algorithm=None):
        """
        Overload C++ constructor for faster init.
        """
        super(PairwiseAligner, self).__init__(scoring, gap_penalty)
        if target is not None:
            self.target = target
        if query is not None:
            self.query = query
        if algorithm is not None:
            self.algorithm = algorithm