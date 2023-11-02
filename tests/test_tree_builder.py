"""
Unit tests for tree builder.
"""
import os
import sys
import tempfile

from sequence_analysis import SeqSet, Sequence, TreeBuilder

sys.path.insert(0, '.')

class TestTreeBuilder:
    def test_init(self):
        sset = SeqSet(file_name="aux_files/dna_ex.fasta")
        with tempfile.TemporaryDirectory() as tmpdirname:
            tb = TreeBuilder(sset, tmpdirname)

    def test_build(self):
        sset = SeqSet()
        s1 = Sequence("TGCAACATGTACTACCCCGAGAGATGGATGGACATGTCTAACTACTACATGGATATGCAAGGACGCTACATGGACAGGTGGGGCCGTTACTGCAAC---------------------------------------")
        s1.name = "first"
        s2 = Sequence("CGTAACACGTGTTACCCCGAGAGATGGATGGACATGTCCAACTACTGCATGGATATGCAGGGACGCTACATGGACAGATCGGGCCGTCATTGCAACCCTTCT---------------------------------")
        s2.name = "second"
        s3 = Sequence("CGTCCCACGTGTTACCCCGAGAGATGGATGGACATGTCCAACTACTGCATGGATATGCAGGGACGCTACATGGACAGATCGGGCCGTCATTGCAACCCTTCT---------------------------------")
        s3.name = "third"
        s4 = Sequence("CGTCCCACGTGTTACCCCGAGAGATGGATGGACAGGGCCAACTACTGCATGGATATGCAGGGACGCTACATGGACAGATCGGGCCGTCATTGCAACCCTTCT---------------------------------")
        s4.name = "fourth"
        sset.add_sequence(s1)
        sset.add_sequence(s2)
        sset.add_sequence(s3)
        sset.add_sequence(s4)
        with tempfile.TemporaryDirectory() as tmpdirname:
            tb = TreeBuilder(sset, tmpdirname)

            if tb.iqtree:
                tb.build()
                assert os.path.exists(f"{tmpdirname}/sequences.fasta.treefile")