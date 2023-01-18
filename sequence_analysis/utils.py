"""
This file holds some basic tools and definitions, mostly for sequence I/O.
"""
import re
from Bio.Seq import Seq
import Bio.Align.substitution_matrices
blosum_50 = Bio.Align.substitution_matrices.load(name="BLOSUM50")

aa_alphabet = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

dna_alphabet = ['G','C','T','A']
rna_alphabet = ['G','C','U','A']

# letters in aa alphabet but not in dna or rna
diff_letters = set(aa_alphabet) - (set(dna_alphabet).union(set(rna_alphabet)))

def check_for_pattern(seq_string, regex):
    if not isinstance(seq_string, str):
        print("ERROR: need to provide either a string or a sequence object")
        return None

    # checks for the existence of pattern described by the regex
    p = re.compile(regex)
    m = p.search(seq_string)
    if m is None:
        return False
    else:
        return True

def query_blosum50(aa1, aa2):
    aa1_ind = blosum_50.alphabet.find(aa1)
    aa2_ind = blosum_50.alphabet.find(aa2)

    return blosum_50[aa1_ind,aa2_ind]