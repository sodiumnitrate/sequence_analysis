"""
This file holds the open reading frame class.
"""
from .sequence_analysis_cpp import OpenReadingFrame as OpenReadingFrame_cpp


class OpenReadingFrame(OpenReadingFrame_cpp):
    """
    This class holds the open reading frame information for a given dna or 
    rna sequence.

    Properties:
    ------------
    std::string rna_sequence;
    std::string parent_sequence;
    std::string protein_sequence;
    int start;
    int stop;
    int strand;
    int frame;
    --------------
    """
    pass
