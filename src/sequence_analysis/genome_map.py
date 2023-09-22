"""
GenomeMap class.
"""
from .sequence_analysis_cpp import GenomeMap as GenomeMap_cpp

class GenomeMap(GenomeMap_cpp):
    """
    Python bindings for GenomeMap.

    attributes:
    ------------
    string chromosome_name
    string sample_name

    methods:
    --------
    get_heatmap(start_idx, end_idx)
    """
    pass
