from sequence_analysis import MSA
from sequence_analysis import SeqSet, Sequence

import os
import sys

sys.path.insert(0, '.') 

class TestMSA:
    def test_msa_setup(self):
        seqs = SeqSet()
        seqs.read_fasta("aux_files/cluster_test_2.fasta")

        msa = MSA(seqs)
        assert msa.sequences is not None
        assert msa.input_name is None
        assert msa.output_name is None
        assert msa.op == 1.53
        assert msa.ep == 0.123
        assert msa.aligned_sequences is None
        assert msa.alignment_info is None
        assert msa.elapsed is None
        assert msa.scratch_folder_name == "tmp"

    def test_msa_align(self):
        seqs = SeqSet()
        seqs.read_fasta("aux_files/cluster_test_2.fasta")

        msa = MSA(seqs)
        if msa.mafft:
            msa.align(clean=True)

            assert msa.alignment_info is not None
            assert msa.aligned_sequences is not None
            assert len(msa.aligned_sequences) == len(msa.sequences)
            assert msa.elapsed is not None

            assert not os.path.exists(msa.scratch_folder_name)

    def test_msa_align_by_codon(self):
        seqs = SeqSet()
        s1 = Sequence("ACGATG")
        s1.name = "s1"
        s2 = Sequence("AGGATC")
        s2.name = "s2"
        s3 = Sequence("GGGATCGCG")
        s3.name = "s3"
        seqs.add_records([s1, s2, s3])
        seqs.type = 'dna'

        msa = MSA(seqs)
        if msa.mafft:
            msa.align_by_codon(clean=True)
            assert msa.alignment_info is not None
            assert msa.aligned_sequences is not None
            assert len(msa.aligned_sequences) == len(msa.sequences)