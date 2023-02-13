from sequence_analysis.seq_set import seq_set


class TestSixFrame:
    def test_six_frame_check(self):
        dna = seq_set(file_name="aux_files/sample_dna.fasta")
        seq = dna.records[0]
        true_sequence = seq.six_frame_check("MDM")
        assert (true_sequence is not None)
