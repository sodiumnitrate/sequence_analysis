"""
This file holds some basic tools and definitions, mostly for sequence I/O.
"""
import re
import string
import random
from statistics import mean
import Bio.Align.substitution_matrices
import numpy as np

# the BLOSUM50 matrix
blosum_50 = Bio.Align.substitution_matrices.load(name="BLOSUM50")

# one-letter alphabet for amino acids
aa_alphabet = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H',
               'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# one-letter alphabet for dna
dna_alphabet = ['G', 'C', 'T', 'A']

# one-letter alphabet for rna
rna_alphabet = ['G', 'C', 'U', 'A']

# letters in aa alphabet but not in dna or rna
diff_letters = set(aa_alphabet) - (set(dna_alphabet).union(set(rna_alphabet)))

# Wimley & White's hydrophobicity scale (kcal/mol, delta_G of transfer from POPC surface to water)
# (Wimley & White, Nature Structural Biology, 1996)
ww = {'A': -0.17, 'R': 0, 'N': -0.42, 'D': -1.23, 'C': 0.24, 'Q': -0.58, 'E': -2.02, 'G': -0.01,
      'H': -0.17, 'I': 0.31, 'L': 0.56, 'K': 0, 'M': 0.23, 'F': 1.13, 'P': -0.45, 'S': -0.13,
      'T': -0.14, 'W': 1.85, 'Y': 0.94, 'V': -0.07}

# amino acid molecular weights
# (https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html)
mw_aa = {'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
         'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
         'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
         'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326}


def kcal_mol_to_kbt(energy_in_kcal_mol, temperature):
    """Function to convert energies from kcal/mol to kBT at a given temperature T."""
    return energy_in_kcal_mol * \
        (4184 / 6.0221409e23) / (1.38064852e-23 * temperature)


def check_for_pattern(seq_string, regex):
    """Function that checks for a regex pattern in a given string."""
    if not isinstance(seq_string, str):
        print("ERROR: need to provide either a string or a sequence object")
        return None

    # checks for the existence of pattern described by the regex
    pattern = re.compile(regex)
    match = pattern.search(seq_string)
    if match is None:
        return False
    return True


def query_blosum50(aa1, aa2):
    """Function that returns score from BLOSUM50 given a pair of amino acid letters."""
    aa1_ind = blosum_50.alphabet.find(aa1)
    aa2_ind = blosum_50.alphabet.find(aa2)

    return blosum_50[aa1_ind, aa2_ind]


def find_spans(full_string, matches):
    """Function that finds the spans of the matching substrings in a given string."""
    chars_absent_from_full_string = list(
        set(string.ascii_uppercase) - set(full_string))
    dummy_char = random.choice(chars_absent_from_full_string)

    spans = []
    for match in matches:
        start_ind = full_string.find(match)
        block = dummy_char * len(match)
        full_string = full_string[0:start_ind] + \
            block + full_string[start_ind + len(match):]
        spans.append((start_ind, start_ind + len(match) - 1))

    return spans


def invert_spans(spans):
    """Function that returns a list spans in between spans in a given list."""
    inverted = []
    for i in range(len(spans) - 1):
        start = spans[i][1] + 1
        end = spans[i + 1][0] - 1
        inverted.append((start, end))

    return inverted


def add_dicts(dict1, dict2):
    """
    Function that adds two dictionaries with identical keys such that
    the values corresponding to the same key get added.
    """
    # adds together values in dicts given the keys are exactly the same
    assert (dict1.keys() == dict2.keys())

    for key in dict1.keys():
        dict1[key] += dict2[key]

    return dict1


def movmean(nums, window=5):
    """Function that calculates the moving mean of a list of numbers."""
    # implement a simple moving mean
    assert (window > 0)
    assert (isinstance(window, int))

    # check if nums is a list of numbers
    if np.isscalar(nums):
        print("ERROR: you must provide a list or array of numbers for calculating moving mean.")
        return None

    # determine radius
    if window % 2 != 0:
        left_r = window // 2
        right_r = window // 2
    else:
        # matlab style
        left_r = window // 2
        right_r = window // 2 - 1

    averaged = []
    for i in range(len(nums)):
        start = max(i - left_r, 0)
        end = min(i + right_r, len(nums))

        averaged.append(mean(nums[start:end + 1]))

    return averaged


def gen_non_overlapping_points(n_points, min_dist, box_length, max_iter=100):
    """Function that generates coordinates
    for points that are separated by at least min_dist."""
    coords = []
    min_dist_2 = min_dist ** 2
    for i in range(n_points):
        placed = False
        for k in range(max_iter):
            x = np.random.uniform(0, box_length)
            y = np.random.uniform(0, box_length)
            if len(coords) == 0:
                coords.append([x, y])
                placed = True
                break

            for coord in coords:
                clash = False
                dist2 = (coord[0] - x) ** 2 + (coord[1] - y) ** 2
                if dist2 < min_dist_2:
                    clash = True
                    break
            if not clash:
                placed = True
                coords.append([x, y])
                break
        if not placed:
            print(
                f"ERROR: node not placed in max_iter={max_iter}. Increase max_iter or decrease min_dist.")
            return None

    return coords


def generate_random_color(n_colors):
    """Function that generates list of random colors."""
    colors = []
    for i in range(n_colors):
        hexadecimal = [
            "#" + ''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
        colors.append(hexadecimal[0])

    return colors
