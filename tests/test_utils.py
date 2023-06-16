import sequence_analysis.utils as sa_utils
from sequence_analysis.sequence import sequence
from Bio.Seq import Seq
import numpy as np

import pdb

class TestUtils:
    def test_blosum_query(self):
        assert (sa_utils.query_blosum50('A', 'A') == 5)
        assert (sa_utils.query_blosum50('C', 'H') == -3)
        assert (sa_utils.query_blosum50('K', 'L') == -3)

    def test_check_for_pattern(self):
        # check string
        s = 'ALALKALIALMMMALKALIMMM'
        assert (sa_utils.check_for_pattern(s, "MMM"))
        assert (sa_utils.check_for_pattern(s, "M+"))

        # check sequence obj
        s = sequence(seq=s, seq_type='protein')
        assert (s.check_for_pattern("MMM"))
        assert (s.check_for_pattern("M+"))

    def test_check_for_pattern_wrong_input(self):
        s = 50
        assert (sa_utils.check_for_pattern(s, "MMM") is None)

    def test_add_dicts(self):
        dict1 = {0: 1, 1: 2}
        dict2 = {0: 1, 1: 5}

        dict3 = sa_utils.add_dicts(dict1, dict2)
        assert (dict3[0] == 2)
        assert (dict3[1] == 7)

    def test_movmean(self):
        nums = [1, 2, 5, 3, 10, 5]
        averaged = sa_utils.movmean(nums, window=2)
        assert (averaged == [1, 1.5, 3.5, 4, 6.5, 7.5])

        averaged = sa_utils.movmean(nums, window=3)
        assert (averaged == [1.5, 8 / 3, 10 / 3, 6, 6, 7.5])

        averaged = sa_utils.movmean(nums, window=20)
        assert (averaged == [13 / 3, 13 / 3, 13 / 3, 13 / 3, 13 / 3, 13 / 3])

        averaged = sa_utils.movmean(5, window=10)
        assert (averaged is None)

    def test_gen_non_overlapping_points(self):
        coords = sa_utils.gen_non_overlapping_points(10, 0.5, 10)
        assert (coords is not None)
        assert (len(coords) == 10)

    def test_sort_array_by_column(self):
        test = [['a',0,5,'aa'],['b',1,0,'bb'],['c',2,1,'cc'],['d',3,5,'dd']]

        sorted_test = sa_utils.sort_array_by_column(test, 1)
        assert len(sorted_test) == len(test)

        sorted_test = sa_utils.sort_array_by_column(test, 0)
        assert sorted_test[0][0] == 'a'

        sorted_test = sa_utils.sort_array_by_column(test, 5)
        assert sorted_test is None

        test = [['a',0], ['b',1,2]]
        sorted_test = sa_utils.sort_array_by_column(test, 0)
        assert sorted_test is None
    
    def test_split_string_into_list(self):
        seq = "A" * 1000
        result = sa_utils.split_string_into_list(seq, 80)

        assert len(result) == 1000 // 80 + 1
        assert len(result[0]) == 80

    def test_check_range(self):
        range_idx = (0,10)
        assert sa_utils.check_range(range_idx)

        range_idx = [0,50]
        assert sa_utils.check_range(range_idx)

        range_idx = np.array([0,50])
        assert sa_utils.check_range(range_idx)

        range_idx = [5]
        assert not sa_utils.check_range(range_idx)

    def test_is_idx_in_range(self):
        idx = 5
        range_idx = (0,10)
        assert sa_utils.is_idx_in_range(idx, range_idx)

        idx = 10
        range_idx = (0,10)
        assert not sa_utils.is_idx_in_range(idx, range_idx)

        idx = 50
        range_idx = (1,10)
        assert not sa_utils.is_idx_in_range(idx, range_idx)

        idx = 0
        range_idx = (2,10)
        assert not sa_utils.is_idx_in_range(idx, range_idx)

    def test_color_substring(self):
        test_string = "LOREM IPSUM DOLOR"
        new_string = sa_utils.color_substring(test_string, (6,11))

        # TODO: is this robust?
        assert new_string[6] == '\x1b'
        assert new_string[16] == "\x1b"

    def test_select_substring_with_gaps(self):
        a = "AAAABBBAAAAAA"
        b = sa_utils.select_substring_with_gaps(a, 4)
        assert b == '------AAAABBBAAAAAA----'

        b = sa_utils.select_substring_with_gaps(a, 4, n_left=2, n_right=2)
        assert b == 'AABBBAA'

    def test_write_in_columns(self):
        a = "A" * 500
        f = open('aux_files/write_in_cols.dat','w')
        sa_utils.write_in_columns(f, a, ncols=10)
        f.close()

        with open('aux_files/write_in_cols.dat', 'r') as f:
            lines = f.readlines()

        assert len(lines) == 50

    def test_write_in_columns_2(self):
        a = "A" * 500
        f = open('aux_files/write_in_cols.dat','w')
        sa_utils.write_in_columns(f, a, ncols=79)
        f.close()

        with open('aux_files/write_in_cols.dat', 'r') as f:
            lines = f.readlines()

        assert len(lines) == 7
        assert len(lines[-1].strip()) == 26

    def test_find_kmers(self):
        test_string = "ACACACA"
        two_mers = sa_utils.find_kmers_in_string(test_string, 2)
        three_mers = sa_utils.find_kmers_in_string(test_string, 3)

        assert len(two_mers) == 2
        assert len(three_mers) == 2
        assert len(list(three_mers.keys())[0]) == 3
        assert len(list(three_mers.keys())[1]) == 3

    def test_find_kmers_2(self):
        test_string = "ALKALI"
        five_mers = sa_utils.find_kmers_in_string(test_string, 5)
        assert len(list(five_mers.keys())) == 2
        assert five_mers['ALKAL'] == 1
        assert five_mers['LKALI'] == 1

        six_mers = sa_utils.find_kmers_in_string(test_string, 6)
        assert len(list(six_mers.keys())) == 1
        assert six_mers['ALKALI'] == 1

        seven = sa_utils.find_kmers_in_string(test_string, 7)
        assert not seven

    def test_find_kmers_in_list(self):
        seq = list("AAACAAACALKAA")
        twomers = sa_utils.find_kmers_in_list(seq, 2)
        assert len(twomers) == 6
        assert twomers[('A','A')] == 5

        threemers = sa_utils.find_kmers_in_list(seq, 3)
        assert len(threemers) == 8
        assert threemers[('A','A','A')] == 2

    def test_check_range_overlap(self):
        range_1 = (0, 10)
        range_2 = (5, 500)

        assert sa_utils.check_range_overlap(range_1, range_2)

        range_1 = (0, 10)
        range_2 = (0, 10)

        assert sa_utils.check_range_overlap(range_1, range_2)

        range_1 = (5, 50)
        range_2 = (55, 100)

        assert not sa_utils.check_range_overlap(range_1, range_2)

        range_1 = 55, 100
        range_2 = 5, 50

        assert not sa_utils.check_range_overlap(range_1, range_2)