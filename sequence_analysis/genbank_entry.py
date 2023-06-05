"""
This file contains the GenBankEntry class.
"""
from Bio import Entrez
from sequence_analysis.sequence import sequence

class GenBankEntry:
    """
    Class that holds info about a genbank entry.
    """
    def __init__(self,
                 accession_code,
                 percent_identity=None,
                 e_value=None,
                 query=None,
                 db='nucleotide'):
        """
        Function to initialize GenBankEntry object.
        """
        self.accession_code = accession_code
        # these can be obtained from blast alignments
        self.percent_identity = percent_identity
        self.e_value = e_value
        self.query = query

        self.db = db

        # attributes to get from genbank
        self.text = None
        self.locus = None
        self.definition = None
        self.organism = None
        self.dna_sequence = None
        self.protein_sequence = None
        self.protein_name = None
        

    def fetch(self, e_mail, skip_origin=True):
        """
        Function to fetch info from genbank.
        """
        Entrez.email = e_mail
        handle = Entrez.efetch(db=self.db,
                               id=self.accession_code,
                               rettype='gb')

        # TODO: some sort of an exception handling here
        text = handle.read()
        self.text = text

        self.get_locus()
        self.get_definition()
        self.get_organism()
        if not skip_origin:
            self.get_origin()
        if self.db == "nucleotide":
            self.get_protein_nucleotide()
        elif self.db == "protein":
            self.get_protein_name_protein()

    def extract_main_keyword(self, keyword):
        """
        Function to extract info given a keyword.
        """
        if self.text is None:
            print("ERROR: run fetch() first.")
            return None

        if keyword not in self.text:
            print(f"WARNING: unable to find {keyword} in data.")
            return None

        return self.text.split(keyword)[1].split('\n')[0].strip()

    def get_locus(self):
        """Function to get LOCUS info."""
        self.locus = self.extract_main_keyword("LOCUS")

    def get_definition(self):
        """Function to get DEFINITION info."""
        self.definition = self.extract_main_keyword("DEFINITION")

    def get_organism(self):
        """Function to get ORGANISM info."""
        organism = self.extract_main_keyword("ORGANISM")
        if organism is None:
            organism = self.extract_main_keyword("SOURCE")

        self.organism = organism

    def get_origin(self):
        """Function to get ORIGIN info."""
        exists = "ORIGIN" in self.text

        if exists:
            string = self.text.split("ORIGIN")[1].split("//")[0].replace("\n","").replace(" ","")
            str2 = "".join(x for x in string if x.isalpha())
            if self.db == "nucleotide":
                self.dna_sequence = sequence(str2)
            elif self.db == "protein":
                self.protein_sequence = sequence(str2)

    def get_protein_nucleotide(self):
        """For a query in the nucleotide db, get protein sequence if it exists."""
        exists = "/translation=" in self.text
        if not exists:
            return

        # make sure that there are not multiple proteins
        n_split = len(self.text.split("/translation="))
        if n_split > 2:
            return

        # handle the case were we have an unnamed protein product
        try:
            string = self.text.split("product=")[1].split("\n")[0].strip('"')
            self.protein_name = string
        except IndexError:
            self.protein_name = None

        string = self.text.split("translation=")[1].split('"')[1].replace("\n","").replace(" ","")
        self.protein_sequence = sequence(string)

    def get_protein_name_protein(self):
        """Function to get protein name if db=protein."""
        exists = "Protein" in self.text
        if not exists:
            return

        string = self.text.split("product=")[1].split("\n")[0].strip('"')
        self.protein_name = string
