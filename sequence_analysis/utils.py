"""
This file holds some basic tools and definitions, mostly for sequence I/O.
"""
import re
from Bio.Seq import Seq
import Bio.Align.substitution_matrices
import string
import random
from statistics import mean
import numpy as np

# the BLOSUM50 matrix
blosum_50 = Bio.Align.substitution_matrices.load(name="BLOSUM50")

# one-letter alphabet for amino acids
aa_alphabet = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

# one-letter alphabet for dna
dna_alphabet = ['G','C','T','A']

# one-letter alphabet for rna
rna_alphabet = ['G','C','U','A']

# letters in aa alphabet but not in dna or rna
diff_letters = set(aa_alphabet) - (set(dna_alphabet).union(set(rna_alphabet)))

# Wimley & White's hydrophobicity scale (kcal/mol, delta_G of transfer from POPC surface to water)
# (Wimley & White, Nature Structural Biology, 1996)
ww = {'A':-0.17, 'R':0, 'N':-0.42, 'D':-1.23, 'C':0.24, 'Q':-0.58, 'E':-2.02, 'G':-0.01, 
'H':-0.17, 'I':0.31,'L':0.56, 'K':0, 'M':0.23, 'F':1.13, 'P':-0.45, 'S':-0.13, 'T':-0.14, 'W':1.85, 'Y':0.94, 'V':-0.07}

# amino acid molecular weights (https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html)
mw_aa = {'A':71.0788, 'R':156.1875, 'N':114.1038, 'D':115.0886, 'C':103.1388,
'E':129.1155, 'Q':128.1307, 'G':57.0519, 'H':137.1411, 'I':113.1594,
'L':113.1594, 'K':128.1741, 'M':131.1926, 'F':147.1766, 'P':97.1167,
'S':87.0782, 'T':101.1051, 'W':186.2132, 'Y':163.1760, 'V':99.1326}

# function to convert energies from kcal/mol to kBT, at a given T
def kcal_mol_to_kbt(en,T):
    return en * (4184/6.0221409e23)/(1.38064852e-23*T)

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

def find_spans(full_string, matches):
    p = list(set(string.ascii_uppercase) - set(full_string))
    c = random.choice(p)

    spans = []
    for match in matches:
        start_ind = full_string.find(match)
        block = c * len(match)
        full_string = full_string[0:start_ind] + block + full_string[start_ind+len(match):]
        spans.append((start_ind,start_ind+len(match)-1))

    return spans

def invert_spans(spans):
    inverted = []
    for i in range(len(spans)-1):
        start = spans[i][1] + 1
        end = spans[i+1][0] - 1
        inverted.append((start,end))

    return inverted

def add_dicts(dict1, dict2):
    # adds together values in dicts given the keys are exactly the same
    assert(dict1.keys() == dict2.keys())

    for key in dict1.keys():
        dict1[key] += dict2[key]

    return dict1

def movmean(nums,window=5):
    # implement a simple moving mean
    assert(window > 0)
    assert(isinstance(window,int))

    # check if nums is a list of numbers
    if np.isscalar(nums):
        print("ERROR: you must provide a list or array of numbers for calculating moving mean.")
        return None

    # determine radius
    if window % 2 != 0:
        lr = window // 2
        rr = window // 2
    else:
        # matlab style
        lr = window // 2
        rr = window // 2 - 1

    averaged = []
    for i in range(len(nums)):
        start = max(i-lr,0)
        end = min(i+rr,len(nums))

        averaged.append(mean(nums[start:end+1]))

    return averaged