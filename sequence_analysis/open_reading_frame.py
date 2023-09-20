"""
This file holds the open reading frame class.
"""
import sequence_analysis_cpp


class OpenReadingFrame(sequence_analysis_cpp.OpenReadingFrame):
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