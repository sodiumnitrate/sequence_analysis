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
    pass
