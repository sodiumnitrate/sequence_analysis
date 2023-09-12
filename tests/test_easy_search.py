"""
Unit tests for easy-search.
"""
import tempfile
from sequence_analysis.easy_search import EasySearch

class TestEasySearch:
    def test_init(self):
        es = EasySearch("aux_files/dna_ex.fasta", "aux_files/dup_test_2.fasta")

        assert es.query_file == "aux_files/dna_ex.fasta"

    def test_set_params(self):
        es = EasySearch("aux_files/dna_ex.fasta", "aux_files/dup_test_2.fasta")
        with tempfile.TemporaryDirectory() as tmpdirname:
            es.set_search_parameters(search_type=3, output_name=f"{tmpdirname}/out.sam", tmp_folder_name=f"{tmpdirname}/tmp")
            assert es.search_type == 3
            assert es.output_name == f"{tmpdirname}/out.sam"
            assert es.cov_mode == 2
            assert es.tmp_folder_name == f"{tmpdirname}/tmp"

    def test_run(self):
        es = EasySearch("aux_files/dna_ex.fasta", "aux_files/dup_test_2.fasta")
        with tempfile.TemporaryDirectory() as tmpdirname:
            es.set_search_parameters(search_type=3, output_name=f"{tmpdirname}/out.sam", tmp_folder_name=f"{tmpdirname}/tmp")

            es.run_search()
            assert es.exec_string is not None