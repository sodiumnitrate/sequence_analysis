"""
Comparing with the results of ORFfinder:
https://www.ncbi.nlm.nih.gov/orffinder/
For human ubiquitin (see aux_files/ubb.fa).
"""
from sequence_analysis.seq_set import seq_set

class TestOrfUbb:
    def test_orf_1(self):
        sequences = seq_set(file_name="aux_files/ubb.fa")
        seq = sequences.records[0]

        orfs = seq.detailed_orf_finder(min_orf_len=70)
        for orf in orfs:
            start = orf.start
            stop = orf.stop
            nuc = orf.rna_sequence.seq
            assert (stop - start) % 3 == 0
            assert len(nuc) % 3 == 0
            prot = orf.protein_sequence.seq
            assert len(prot.strip('*').split('*')) == 1

    def test_orf_2(self):
        sequences = seq_set(file_name='aux_files/ubb.fa')
        prots = seq_set(file_name='aux_files/ubb_prot_orfs.fa')
        seq = sequences.records[0]
        prot_seq = prots.records[0]

        orfs = seq.detailed_orf_finder(min_orf_len=87)

        assert orfs[0].protein_sequence.seq.strip('*') == prot_seq.seq
        assert orfs[0].start == 264
        assert orfs[0].stop == 954
        assert orfs[0].strand == 1
        assert orfs[0].frame == 0 # frame=0 in our code is frame=1 on ncbi

    def test_orf_3(self):
        sequences = seq_set(file_name='aux_files/ubb.fa')
        seq = sequences.records[0]
        orfs = seq.detailed_orf_finder(min_orf_len=87)
        assert len(orfs) == 6
