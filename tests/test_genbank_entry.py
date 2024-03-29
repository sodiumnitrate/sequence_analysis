from sequence_analysis import GenBankEntry
from sequence_analysis import Sequence

class TestGenBankEntry:
    def test_fetch(self):
        accession_code = "HE687196.1"
        entry = GenBankEntry(accession_code=accession_code)
        entry.fetch(e_mail="irem.altan@yale.edu",
                    skip_origin=False)

        assert entry.text is not None
        assert accession_code[:-2] in entry.locus
        assert isinstance(entry.definition, str)
        assert "Sepia officinalis" in entry.organism

        assert isinstance(entry.dna_sequence, Sequence)
        assert len(entry.dna_sequence) > 0

        assert entry.protein_name is not None
        assert isinstance(entry.protein_sequence, Sequence)
        assert len(entry.protein_sequence) > 0

    def test_fetch_2(self):
        accession_code = "AY557479.1"

        entry = GenBankEntry(accession_code=accession_code)
        entry.fetch(e_mail="irem.altan@yale.edu",
                    skip_origin=False)

        assert entry.text is not None
        assert accession_code[:-2] in entry.locus
        assert isinstance(entry.definition, str)
        assert "Doryteuthis pealeii" in entry.organism

        assert isinstance(entry.dna_sequence, Sequence)
        assert len(entry.dna_sequence) > 0

        assert entry.protein_name is None
        assert entry.protein_sequence is None

    def test_fetch_3(self):
        accession_code = "OV121133.1"
        entry = GenBankEntry(accession_code=accession_code)
        entry.fetch(e_mail="irem.altan@yale.edu")

        assert entry.text is not None
        assert accession_code[:-2] in entry.locus
        assert isinstance(entry.definition, str)
        assert "Brassicogethes aeneus" in entry.organism

        assert entry.dna_sequence is None
        assert entry.protein_name is None
        assert entry.protein_sequence is None

    def test_fetch_protein(self):
        code = "CCG28047.1"
        entry = GenBankEntry(accession_code=code, db='protein')
        entry.fetch(e_mail="irem.altan@yale.edu", skip_origin=False)

        assert isinstance(entry.protein_sequence, Sequence)
        assert entry.dna_sequence is None
        assert entry.protein_name is not None

    def test_fetch_protein_2(self):
        code = "BAH14572.1"
        entry = GenBankEntry(accession_code=code, db='protein')
        entry.fetch(e_mail="irem.altan@yale.edu", skip_origin=False)

        assert entry.protein_name is not None

    def test_fetch_protein_3(self):
        code = "5OF9_A"
        entry = GenBankEntry(accession_code=code, db='protein')

        entry.fetch(e_mail="irem.altan@yale.edu", skip_origin=False)

        assert entry.protein_name is not None
