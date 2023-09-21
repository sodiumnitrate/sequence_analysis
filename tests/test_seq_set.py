"""
Unit tests for SeqSet.
"""
import tempfile

from sequence_analysis.seq_set import SeqSet
from sequence_analysis.sequence import Sequence

class TestSeqSet:
    def test_init(self, capsys):
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

        captured = capsys.readouterr()
        # do stuff

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