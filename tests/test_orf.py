from sequence_analysis.sequence import sequence
from sequence_analysis.open_reading_frame import OpenReadingFrame

class TestOrf:
    def test_find_fragments(self):
        seq = sequence('AAACCCATGTTTAAATTTAAATAGAAAGGG')
        orfs = seq.find_fragments(seq,frame=0,order=1,min_orf_len=0)
        assert len(orfs) == 1
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop
        assert parent_str[start:stop] == nuc_str

    def test_find_fragments_2(self):
        seq = sequence('AAACCCATGTTTAAATTTAAATAGAAAGGGA')
        orfs = seq.find_fragments(seq,frame=0,order=1,min_orf_len=0)
        assert len(orfs) == 1
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop
        assert parent_str[start:stop] == nuc_str

    def test_find_fragments_3(self):
        seq = sequence('AAACCCATGTTTAAATTTAAATAGAAAGGGAT')
        orfs = seq.find_fragments(seq,frame=0,order=1,min_orf_len=0)
        assert len(orfs) == 1
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop
        assert parent_str[start:stop] == nuc_str

    def test_find_fragments_4(self):
        seq = sequence('GAAACCCATGTTTAAATTTAAATAGAAAGGGAT')
        seq_fs = seq.frame_shift(frame=1)
        orfs = seq_fs.find_fragments(seq,frame=1,order=1,min_orf_len=0)
        assert len(orfs) == 1
        orfs[0].parent_sequence = seq
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop
        assert parent_str[start:stop] == nuc_str

    def test_find_fragments_5(self):
        seq = sequence('TTTCTATTTAAATTTAAACATGGGTTT')
        seq_fs = seq.reverse_complement().frame_shift(frame=0)
        orfs = seq_fs.find_fragments(seq,frame=0, order=-1, min_orf_len=0)
        assert len(orfs) == 1
        orfs[0].parent_sequence = seq
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop

        comp_str = sequence(parent_str[start:stop], seq_type=seq.type).reverse_complement().seq
        assert comp_str == nuc_str

    def test_find_fragments_6(self):
        seq = sequence('TTTCTATTTAAATTTAAACATGGGTTTA')
        seq_fs = seq.reverse_complement().frame_shift(frame=1)
        orfs = seq_fs.find_fragments(seq,frame=1, order=-1, min_orf_len=0)
        assert len(orfs) == 1
        orfs[0].parent_sequence = seq
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop

        comp_str = sequence(parent_str[start:stop], seq_type=seq.type).reverse_complement().seq
        assert comp_str == nuc_str

    def test_find_fragments_7(self):
        seq = sequence('TTTCTATTTAAATTTAAACATGGGTTTAA')
        seq_fs = seq.reverse_complement().frame_shift(frame=2)
        orfs = seq_fs.find_fragments(seq,frame=2, order=-1, min_orf_len=0)
        assert len(orfs) == 1
        orfs[0].parent_sequence = seq
        nuc_str = orfs[0].rna_sequence.seq
        parent_str = orfs[0].parent_sequence.seq
        start = orfs[0].start
        stop = orfs[0].stop

        comp_str = sequence(parent_str[start:stop], seq_type=seq.type).reverse_complement().seq
        assert comp_str == nuc_str


    def test_check_indices(self):
        seq = sequence('CGCTACGTCTTACGCTGGAGCTCTCATGGATCGGTTCGGTAGGGCTCGATCACATCGCTAGCCAT')
        orfs = seq.detailed_orf_finder(min_orf_len=0)
        for orf in orfs:
            assert (orf.stop - orf.start) % 3 == 0
            prot = orf.protein_sequence.seq
            start = orf.start
            stop = orf.stop

            assert stop - start == len(orf.rna_sequence)
            
            if orf.strand == -1:
                prot_from_ind = sequence(orf.parent_sequence.seq[start:stop],
                                         seq_type=orf.parent_sequence.type).reverse_complement().translate().seq
            else:
                prot_from_ind = sequence(orf.parent_sequence.seq[start:stop],
                                         seq_type=orf.parent_sequence.type).translate().seq
            assert prot == prot_from_ind

    def test_check_indices_with_gaps(self):
        seq = sequence('CGCTACGTCTTACGCTGGAGCTCTCA----TGGATCGGTTCGGTAGGGCTCGATCACATCGCTAGCCAT')
        orfs = seq.detailed_orf_finder(min_orf_len=0)
        for orf in orfs:
            assert (orf.stop - orf.start) % 3 == 0
            prot = orf.protein_sequence.seq
            start = orf.start
            stop = orf.stop

            assert stop - start == len(orf.rna_sequence)
            
            if orf.strand == -1:
                prot_from_ind = sequence(orf.parent_sequence.seq[start:stop],
                                         seq_type=orf.parent_sequence.type).reverse_complement().translate().seq
            else:
                prot_from_ind = sequence(orf.parent_sequence.seq[start:stop],
                                         seq_type=orf.parent_sequence.type).translate().seq
            assert prot == prot_from_ind