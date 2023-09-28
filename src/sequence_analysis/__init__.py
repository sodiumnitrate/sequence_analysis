"""
__init__ file for the sequence_analysis module
"""
from __future__ import annotations
from .sequence_analysis_cpp import __doc__, __version__

from .seq_set import SeqSet
from .sequence import Sequence
from .open_reading_frame import OpenReadingFrame
from .genome_map import GenomeMap
from .genbank_entry import GenBankEntry
from .msa import MSA
from .pairwise_aligner import PairwiseAligner
from .sam_file import SamFile
from .star_map import StarMap
from .tree import Tree
from .easy_search import EasySearch
from .fasta_iterator import FastaIterator

# from .sequence_analysis_cpp import __doc__, __version__, GenomeMap, OpenReadingFrame, PairwiseAligner, SamFile, SeqSet, Sequence

__all__ = ["__doc__", "__version__", "SeqSet", "Sequence", "OpenReadingFrame", "PairwiseAligner", "SamFile", "StarMap", "Tree", "MSA", "GenBankEntry", "GenomeMap", "EasySearch"]
