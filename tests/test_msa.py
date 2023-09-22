from sequence_analysis import MSA
from sequence_analysis import SeqSet

import os

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
        msa.align(clean=True)

        assert msa.alignment_info is not None
        assert msa.aligned_sequences is not None
        assert len(msa.aligned_sequences) == len(msa.sequences)
        assert msa.elapsed is not None

        assert not os.path.exists(msa.scratch_folder_name)
