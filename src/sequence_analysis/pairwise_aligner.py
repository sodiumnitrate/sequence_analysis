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
    def __init__(self, query=None, target=None, algorithm=None):
        """
        Overload C++ constructor for faster init.

        TODO: reconcile this with the overloaded constructor on the C++ side.
        """
        super(PairwiseAligner, self).__init__()
        if query is not None:
            self.query = query
        if target is not None:
            self.target = target
        if algorithm is not None:
            self.algorithm = algorithm