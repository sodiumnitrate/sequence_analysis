"""
Unit tests for SeqSet.
"""
import tempfile

from sequence_analysis import SeqSet
from sequence_analysis import Sequence

import sys
sys.path.insert(0, '.') 

class TestSeqSet:
    def test_init(self):
        sset = SeqSet()
        s1 = Sequence("ACGT")
        s2 = Sequence("GGCCTT")
        sset.records = [s1, s2]
        assert len(sset) == 2
        assert sset.type == 'dna'

        s3 = Sequence("MDMDMDMDYFF")
        sset.records = [s1, s2, s3]
        assert len(sset) == 3
        assert sset.type == 'dna'
        sset.records = [s3, s1, s2]
        assert sset.type == 'protein'
    
    def test_init_2(self):
        sset = SeqSet(list_of_seqs=[Sequence("ACGT"), Sequence("GGCCTT")])
        assert len(sset) == 2
        assert sset.type == "dna"

    def test_init_3(self):
        sset = SeqSet(file_name="aux_files/test.fasta")
        assert len(sset) == 3

    def test_add_sequence(self):
        sset = SeqSet()
        s1 = Sequence("ACGT")
        s2 = Sequence("GGCCTT")
        s3 = Sequence("AAAAAAAA")
        sset.records = [s1]
        assert len(sset) == 1
        sset.add_sequence(s2)
        assert len(sset) == 2
        sset.add_sequence(s3)
        assert len(sset) == 3

    def test_dealign(self):
        sset = SeqSet()
        s1 = Sequence("AC--G-T-")
        s2 = Sequence("GGCCTT--")
        s3 = Sequence("AAAAAAAA")
        sset.records = [s1, s2, s3]
        sset.dealign()
        assert sset.records[0].seq_str == "ACGT"

    def test_alphabetize(self):
        sset = SeqSet()
        s1 = Sequence("ACGTTTT")
        s2 = Sequence("GGCCTTT")
        s3 = Sequence("CCCGACT")
        sset.records = [s1, s2, s3]
        sset.alphabetize()
        assert sset.records[0].seq_str == "ACGTTTT"
        assert sset.records[1].seq_str == s3.seq_str

    def test_remove_duplicates(self):
        sset = SeqSet()
        s1 = Sequence("ACGTTTT")
        s1.name = "s1"
        s2 = Sequence("GGCCTTT")
        s2.name = "s2"
        s3 = Sequence("CCCGACT")
        s3.name = "s3"
        s4 = Sequence("CCCGACT")
        s4.name = "s4"
        sset.records = [s1, s2, s3, s4]
        sset.remove_duplicates()
        assert len(sset) == 3
        assert sset.records[1].name == "s3_s4"

    def test_add_set(self):
        sset = SeqSet()
        sset.records =[Sequence("AAAA"), Sequence("CCCC")]
        assert len(sset) == 2
        sset2 = SeqSet()
        sset2.records = [Sequence("GGGG"), Sequence("TTTT")]
        sset.add_set(sset2)
        assert len(sset) == 4

    def test_write_fasta(self):
        sset = SeqSet()
        sset.records = [Sequence("GGGG"), Sequence("TTTT"), Sequence("C"*100)]
        with tempfile.TemporaryDirectory() as tmpdirname:
            sset.write_fasta(f"{tmpdirname}/out.fasta")
            ct = 0
            ln_ct = 0
            with open(f"{tmpdirname}/out.fasta",'r') as f:
                for line in f:
                    if line.startswith(">"):
                        ct += 1
                    ln_ct += 1

        assert ct == 3
        assert ln_ct == 7

    def test_read_fasta(self):
        sset = SeqSet()
        sset.read_fasta("aux_files/test.fasta")
        assert len(sset) == 3
        assert sset.records[0].name == 'seq0'
        assert sset.records[1].name == "seq1"
        assert sset.records[2].name == ""
        assert sset.records[2].seq_str == "AAAAA"

    def test_read_fasta_2(self):
        sset = SeqSet()
        sset.read_fasta("aux_files/reference_set_dna.fasta")
        assert len(sset) == 3
        assert sset.records[-1].name == "UBIQUITIN_ESC_2"

    def test_write_phylip(self):
        sset = SeqSet()
        s1 = Sequence("A"*230)
        s1.name = 's1'
        s2 = Sequence("C"*230)
        s2.name = 's2'
        s3 = Sequence("T"*230)
        s3.name = 's3'
        sset.records = [s1, s2, s3]
        assert sset.records[0].name == 's1'

        sset.write_phylip("aux_files/test_s.phy", "sequential")
        sset.write_phylip("aux_files/test_i.phy", "interleaved")

    def test_slice(self):
        sset = SeqSet(file_name="aux_files/test.fasta")
        first = sset[0]
        assert isinstance(first.seq_str, str)
        last_two = sset[1:]
        assert len(last_two) == 2

    def test_trim_gaps(self):
        s1 = Sequence("--ACGT---")
        s2 = Sequence("-AACGT---")
        s3 = Sequence("-AACGTA--")
        sset = SeqSet(list_of_seqs=[s1, s2, s3])
        sset.trim_gaps()

        for seq in sset:
            assert '-' not in seq.seq_str
            assert len(seq) == 4

    def test_pairwise_distance(self):
        s1 = Sequence("ACGTACGT")
        s2 = Sequence("GGGGGG")
        s3 = Sequence("CGCCCCCC")
        s4 = Sequence("CCCCCCC")

        sset = SeqSet(list_of_seqs=[s1, s2, s3, s4])
        adjacency = sset.pairwise_distance()
        for i in range(len(adjacency)):
            assert adjacency[i][i] / len(sset[i]) == 1
            for j in range(len(adjacency)):
                assert adjacency[i][j] == adjacency[j][i]

    def test_find_subset_with_names(self):
        s1 = Sequence("ACGTACGT")
        s1.name = "s1"
        s2 = Sequence("GGGGGG")
        s2.name = "s2"
        s3 = Sequence("CGCCCCCC")
        s3.name = "s3"
        s4 = Sequence("CCCCCCC")
        s4.name = "s4"

        sset = SeqSet(list_of_seqs=[s1, s2, s3, s4])
        newset = sset.find_subset_with_names(["s2", "s1"])
        assert len(newset) == 2
        names = [s.name for s in newset.records]
        assert "s1" in names
        assert "s2" in names

    def test_trim_gaps(self):
        s1 = Sequence("ACGTACGT")
        s1.name = "s1"
        s2 = Sequence("GGGGGG--")
        s2.name = "s2"
        s3 = Sequence("CGCCCCCC")
        s3.name = "s3"
        s4 = Sequence("CCCCCCC-")
        s4.name = "s4"

        sset = SeqSet(list_of_seqs=[s1, s2, s3, s4])

        sset.trim_gaps()
        assert len(sset) == 4
        for seq in sset:
            assert len(seq) == 6

        sset_2 = SeqSet(list_of_seqs=[s1, s2, s3, s4])
        sset_2.trim_gaps(0.7)
        assert len(sset_2) == 4
        for seq in sset_2:
            assert len(seq) == 7

    def test_get_seq_with_name(self):
        s1 = Sequence("ACGTACGT")
        s1.name = "s1"
        s2 = Sequence("GGGGGG--")
        s2.name = "s2"
        s3 = Sequence("CGCCCCCC")
        s3.name = "s3"
        s4 = Sequence("CCCCCCC-")
        s4.name = "s4"

        sset = SeqSet(list_of_seqs=[s1, s2, s3, s4])
        selected = sset.get_sequence_with_name('s4')
        assert selected.name == 's4'