import sequence_analysis.utils as sa_utils
from sequence_analysis.sequence import sequence
from Bio.Seq import Seq

class TestUtils:
    def test_blosum_query(self):
        assert(sa_utils.query_blosum50('A','A') == 5)
        assert(sa_utils.query_blosum50('C','H') == -3)
        assert(sa_utils.query_blosum50('K','L') == -3)

    def test_check_for_pattern(self):
        # check string
        s = 'ALALKALIALMMMALKALIMMM'
        assert( sa_utils.check_for_pattern(s, "MMM"))
        assert( sa_utils.check_for_pattern(s, "M+"))

        # check sequence obj
        s = sequence(seq=s, type='protein')
        assert( sa_utils.check_for_pattern(s, "MMM"))
        assert(sa_utils.check_for_pattern(s, "M+"))

        # check biopythons seq object
        s = Seq(s.seq)
        assert( sa_utils.check_for_pattern(s, "MMM"))
        assert(sa_utils.check_for_pattern(s, "M+"))
