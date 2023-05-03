from sequence_analysis.genbank_entry import GenBankEntry
from sequence_analysis.sequence import sequence

import pdb

class TestGenBankEntry:
    def test_fetch(self):
        accession_code = "HE687196.1"
        entry = GenBankEntry(accession_code=accession_code)
        entry.fetch(e_mail="irem.altan@yale.edu")

        assert entry.text is not None
        assert accession_code[:-2] in entry.locus
        assert isinstance(entry.definition, str)
        assert "Sepia officinalis" in entry.organism

        assert isinstance(entry.dna_sequence, sequence)
        assert len(entry.dna_sequence) > 0

        assert entry.protein_name is not None
        assert isinstance(entry.protein_sequence, sequence)
        assert len(entry.protein_sequence) > 0