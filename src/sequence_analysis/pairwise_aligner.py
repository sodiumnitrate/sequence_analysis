"""
This file holds the PairwiseAligner class and related methods.
"""
from .sequence_analysis_cpp import PairwiseAligner as PairwiseAligner_cpp
from math import ceil

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

    def print(self, char_per_line=80):
        """
        Function to print alignment.
        """
        query_aligned = self.get_query_aligned()
        target_aligned = self.get_target_aligned()
        match = self.get_match_string()
        n_lines = int(ceil(len(match) / char_per_line))
        for i in range(n_lines):
            print(query_aligned[i*char_per_line:(i+1)*char_per_line])
            print(match[i*char_per_line:(i+1)*char_per_line])
            print(target_aligned[i*char_per_line:(i+1)*char_per_line])


    # TODO: refactor below.
    def get_query_range(self):
        """
        Function to get the range of the query that aligns with the target.
        """
        query_aligned = self.get_query_aligned().replace('-','')
        start = self.query.find(query_aligned)
        end = self.query[::-1].find(query_aligned[::-1])
        end = len(self.query) - end - 1
        return start, end

    def get_target_range(self):
        """
        Function to get the range of the target that aligns with the query.
        """
        target_aligned = self.get_target_aligned().replace('-','')
        start= self.target.find(target_aligned)
        end = self.target[::-1].find(target_aligned[::-1])
        end = len(self.target) - end - 1
        return start, end
