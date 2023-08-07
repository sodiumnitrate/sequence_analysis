"""
Unit tests for the seq_match object.
"""

from sequence_analysis.seq_match import SeqMatch, Match

class TestSeqMatch:
    def test_init_match(self):
        test_dict = {"query":"query", "target":"target", "evalue":0.0}
        m = Match(target="target",query="query",properties=test_dict)

        assert m.target == "target"

    def test_init_seq_match(self):
        files = ["aux_files/match/part-00070", "aux_files/match/part-00462"]
        m = SeqMatch(files)

        assert m.strip_chars == "()'"

    def test_seq_match_read(self):
        files = ["aux_files/match/part-00070", "aux_files/match/part-00462"]
        m = SeqMatch(files)
        m.read_data()

        assert len(m.matches) == 20

    def test_seq_match_read_2(self):
        file = "aux_files/match/part-00070"
        m = SeqMatch(file)
        m.read_data()

        assert len(m.matches) == 15
        assert isinstance(m.matches[0], Match)

    def test_index_data(self):
        file = "aux_files/match/part-00070"
        m = SeqMatch(file)
        m.read_data()

        assert m.query_index is not None
        assert m.target_index is not None

        assert sum([len(val) for val in m.target_index.values()]) == len(m.matches)
        assert sum([len(val) for val in m.query_index.values()]) == len(m.matches)

    def test_load_sequences(self):
        file = "aux_files/match/test"
        m = SeqMatch(file, keywords=["query","target","evalue"])
        m.read_data()

        m.load_sequences("aux_files/cluster_test.fasta","aux_files/cluster_test_2.fasta")

        assert m.query_sequences is not None
        assert m.target_sequences is not None
        