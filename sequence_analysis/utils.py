"""
This file holds some basic tools and definitions, mostly for sequence I/O.
"""

aa_alphabet = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

dna_alphabet = ['G','C','T','A']
rna_alphabet = ['G','C','U','A']

# letters in aa alphabet but not in dna or rna
diff_letters = set(aa_alphabet) - (set(dna_alphabet).union(set(rna_alphabet)))