from sequence_analysis.genbank_entry import GenBankEntry
from sequence_analysis.sequence import sequence

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

        assert isinstance(entry.dna_sequence, sequence)
        assert len(entry.dna_sequence) > 0

        assert entry.protein_name is not None
        assert isinstance(entry.protein_sequence, sequence)
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

        assert isinstance(entry.dna_sequence, sequence)
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