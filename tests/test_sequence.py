"""
Unit tests for sequence.

TODO:
- more edge cases for type assignment
- more translate tests
"""
import tempfile
from sequence_analysis.sequence import Sequence

class TestSequence:
    def test_init(self):
        seq = Sequence("ACGT")
        assert seq.seq_str == "ACGT"
        seq.name = "test"
        assert seq.name == "test"
        seq.type = "dna"
        assert seq.type == "dna"

        seq = Sequence("MQIFYQWMDPP")
        assert seq.type == "protein"

    def test_slice(self):
        seq = Sequence("ACGTACGT")
        sliced = seq[1:4]
        assert sliced == "CGT"

    def test_set_type(self):
        seq = Sequence("ACGTCCGGAATT")
        seq.set_type()
        assert seq.type == "dna"

        seq.set_type("blah")
        assert seq.type == "blah"

    def test_find_codon(self):
        seq = Sequence("ACGTCCGGAATTG")
        pos = seq.find_codon("ATT")
        assert pos == 9

        pos = seq.find_codon("CGT")
        assert pos == -1

    def test_reverse(self):
        seq = Sequence("ACGT")
        seq_rev = seq.reverse()
        assert seq_rev.seq_str == seq.seq_str[::-1]

    def test_complement(self):
        seq = Sequence("ACGT")
        assert seq.seq_str == "ACGT"
        seq_comp = seq.complement()
        assert seq_comp.seq_str == "TGCA"

    def test_reverse_complement(self):
        seq = Sequence("ACGTT")
        seq_rc = seq.reverse_complement()
        assert seq_rc.seq_str == "AACGT"

    def test_frame_shift(self):
        seq = Sequence("ACGTACGTT")
        shifted = seq.frame_shift(0)
        assert seq.seq_str == shifted.seq_str
        shifted = seq.frame_shift(1)
        assert seq.seq_str[1:] == shifted.seq_str
        shifted = seq.frame_shift(2)
        assert seq.seq_str[2:] == shifted.seq_str

    def test_translate(self):
        seq = Sequence("atggcggatcagctgaccgaagaacagattgcggaatttaaagaagcgtttagcctgtttgataaagatggcgatggcaccattaccaccaaagaactgggcaccgtgatgcgcagcctgggccagaacccgaccgaagcggaactgcaggatatgattaacgaagtggatgcggatggcaacggcaccattgattttccggaatttctgaccatgatggcgcgcaaaatgaaagataccgatagcgaagaagaaattcgcgaagcgtttcgcgtgtttgataaagatggcaacggctatattagcgcggcggaactgcgccatgtgatgaccaacctgggcgaaaaactgaccgatgaagaagtggatgaaatgattcgcgaagcggatattgatggcgatggccaggtgaactatgaagaatttgtgcagatgatgaccgcgaaa")
        calmodulin = seq.translate()
        assert calmodulin.type == 'protein'
        assert seq.type == 'dna'

        assert calmodulin.seq_str == "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK"

    def test_write_fasta(self):
        seq = Sequence("MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK")
        seq.name = "calmodulin"
        with tempfile.TemporaryDirectory() as tmpdirname:
            seq.write_fasta(f"{tmpdirname}/out.fasta")
            with open(f"{tmpdirname}/out.fasta", 'r') as f:
                first = f.readline()
                second = f.readline()
                assert first.startswith(">")
                assert seq.name in first
                assert len(second.strip()) == 79
                assert second.startswith("MADQL")

    def test_get_orf(self):
        prot_seq = Sequence("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGC")
        assert prot_seq.type == "protein"
        ubiquitin = Sequence("AGTGCCGGCGCTGCAAGGAAGTTTCCAGAGCTTTCGAGGAAGGTTTCTTCAACTCAAATTCATCCGCCTGATAATTTTCTTATATTTTCCTAAAGAAGGAAGAGAAGCGCATAGAGGAGAAGGGAAATAATTTTTTAGGAGCCTTTCTTACGGCTATGAGGAATTTGGGGCTCAGTTGAAAAGCCTAAACTGCCTCTCGGGAGGTTGGGCGCGGCGAACTACTTTCAGCGGCGCACGGAGACGGCGTCTACGTGAGGGGTCAAAATGCAGATCTTCGTGAAAACCCTTACCGGCAAGACCATCACCCTTGAGGTGGAGCCCAGTGACACCATCGAAAATGTGAAGGCCAAGATCCAGGATAAGGAAGGCATTCCCCCCGACCAGCAGAGGCTCATCTTTGCAGGCAAGCAGCTGGAAGATGGCCGTACTCTTTCTGACTACAACATCCAGAAGGAGTCGACCCTGCACCTGGTCCTGCGTCTGAGAGGTGGTATGCAGATCTTCGTGAAGACCCTGACCGGCAAGACCATCACCCTGGAAGTGGAGCCCAGTGACACCATCGAAAATGTGAAGGCCAAGATCCAGGATAAAGAAGGCATCCCTCCCGACCAGCAGAGGCTCATCTTTGCAGGCAAGCAGCTGGAAGATGGCCGCACTCTTTCTGACTACAACATCCAGAAGGAGTCGACCCTGCACCTGGTCCTGCGTCTGAGAGGTGGTATGCAGATCTTCGTGAAGACCCTGACCGGCAAGACCATCACTCTGGAGGTGGAGCCCAGTGACACCATCGAAAATGTGAAGGCCAAGATCCAAGATAAAGAAGGCATCCCCCCCGACCAGCAGAGGCTCATCTTTGCAGGCAAGCAGCTGGAAGATGGCCGCACTCTTTCTGACTACAACATCCAGAAAGAGTCGACCCTGCACCTGGTCCTGCGCCTGAGGGGTGGCTGTTAATTCTTCAGTCATGGCATTCGCAGTGCCCAGTGATGGCATTACTCTGCACTATAGCCATTTGCCCCAACTTAAGTTTAGAAATTACAAGTTTCAGTAATAGCTGAACCTGTTCAAAATGTTAATAAAGGTTTCGTTGCATGGTA")
        ubiquitin.name = "ubiquitin"
        assert ubiquitin.type == "dna"
        orfs = ubiquitin.get_open_reading_frames(min_len=87)
        assert len(orfs) == 6

        expasy = ["MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGC",
        "MAFAVPSDGITLHYSHLPQLKFRNYKFQ",
        "MSLCWSGGMPSLSWILAFTFSMVSLGSTSRVMVLPVRVFTKICIPPLRRRTRCRVDSFWML",
        "MSLCWSGGMPSLSWILAFTFSMVSLGSTSRVMVLPVRVFTKICIPPLRRRTRCRVDSFWML",
        "MSLCWSGGMPSLSWILAFTFSMVSLGSTSRVMVLPVRVFTKICILTPHVDAVSVRR",
        "MPSLGTANAMTEELTATPQAQDQVQGRLFLDVVVRKSAAIFQLLACKDEPLLVGGDAFFILDLGLHIFDGVTGLHLQSDGLAGQGLHEDLHTTSQTQDQVQGRLLLDVVVRKSAAIFQLLACKDEPLLVGRDAFFILDLGLHIFDGVTGLHFQGDGLAGQGLHEDLHTTSQTQDQVQGRLLLDVVVRKSTAIFQLLACKDEPLLVGGNAFLILDLGLHIFDGVTGLHLKGDGLAGKGFHEDLHFDPSRRRRLRAPLKVVRRAQPPERQFRLFN"]

        ours = [orf.protein_sequence.strip("*") for orf in orfs]

        assert set(expasy) == set(ours)