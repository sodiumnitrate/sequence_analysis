import sequence_analysis.utils as sa_utils
from sequence_analysis.sequence import sequence
from Bio.Seq import Seq


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
