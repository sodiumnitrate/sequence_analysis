"""
File that defines a class for organizing the result of mmseqs search or blast.

#TODO: more robust file parsing. XML support.
#TODO: parallelize?

See pg. 54 for the keywords in https://mmseqs.com/latest/userguide.pdf
"""

import copy

from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence
from sequence_analysis.pairwise_alignment import pairwise_alignment
from sequence_analysis.utils import file_name_check


default_keywords = ["query", "target", "evalue", "qlen", "qstart", "qend",
                    "tstart", "tend", "qcov", "raw"]
standard_blast_keywords = ["query", "target", "evalue", "gapopen",
                           "pident", "fident", "nident", "qstart",
                           "qend", "qlen", "tstart", "tend", "tlen",
                           "alnlen", "raw", "bits", "cigar", "qseq",
                           "tseq", "qaln", "taln", "qheader", "theader",
                           "qframe", "tframe", "mismatch", "qcov",
                           "tcov", "empty", "taxid", "taxname",
                           "taxlineage", "qset", "qsetid", "tset", "tsetid"]
                           
standard_blast_keyword_types = {"query":str, "target":str, "evalue":float, "gapopen":int,
                           "pident":float, "fident":float, "nident":int, "qstart":int,
                           "qend":int, "qlen":int, "tstart":int, "tend":int, "tlen":int,
                           "alnlen":int, "raw":float, "bits":float, "cigar":str, "qseq":str,
                           "tseq":str, "qaln":str, "taln":str, "qheader":str, "theader":str,
                           "qframe":int, "tframe":int, "mismatch":int, "qcov":float,
                           "tcov":float, "empty":str, "taxid":str, "taxname":str,
                           "taxlineage":str, "qset":str, "qsetid":int, "tset":str, "tsetid":int}

class Match:
    def __init__(self, target, query, properties):
        self.target = target
        self.query = query
        self.properties = properties

    def check_data(self):
        """Check the validity of input data."""
        if not isinstance(self.target, str):
            raise TypeError
        if not isinstance(self.query, str):
            raise TypeError

        if not isinstance(self.properties, dict):
            raise TypeError

        for key, val in self.properties.items():
            if key not in standard_blast_keywords:
                raise ValueError
            if not isinstance(val, standard_blast_keyword_types[key]):
                raise TypeError


class SeqMatch:
    def __init__(self, file_names, keywords=default_keywords, delimiter=",", strip_chars="()'"):
        self.file_names = file_names
        self.keywords = keywords
        self.delimiter = delimiter
        self.strip_chars = strip_chars

        self.check_input_values()

        self.query_names = None
        self.target_names = None
        self.matches = None

        self.query_index = None
        self.target_index = None

        self.query_sequences = None
        self.target_sequences = None

        self.target_truncated = None

    def check_input_values(self):
        """
        Function to check whether the values the user provided make sense.
        """
        # file_name can be a string for a single file, or a list of strings for multiple files
        if not isinstance(self.file_names, str):
            if not isinstance(self.file_names, list):
                print(f"ERROR: file_names of {type(self.file_names)} not supported.")
                raise TypeError
            else:
                if not all([isinstance(f, str) for f in self.file_names]):
                    print(f"ERROR: file_names is given as a list, but the list doesn't contain all strings.")
                    raise TypeError

        # keywords should be provided as a list, and should match the standard blast keywords
        # TODO: add support for a single string?
        if not isinstance(self.keywords, list):
            raise TypeError

        if not all([isinstance(f, str) for f in self.keywords]):
            raise TypeError

        for keyword in self.keywords:
            if keyword not in standard_blast_keywords:
                print(f"ERROR: {keyword} not a standard keyword.")
                raise ValueError

        # query and target must be in the keywords (otherwise the match is meaningless)
        if "query" not in self.keywords or "target" not in self.keywords:
            print("ERROR: query and target must be in the keywords.")
            raise ValueError

        # delimiter must be a single character
        if len(self.delimiter) != 1:
            print(f"ERROR: {self.delimiter} is not a valid delimiter. Use a single character.")
            raise ValueError

    def check_files(self):
        """
        Function to check that we can open the listed files and they
        indeed have the correct number of fields.
        """
        if not isinstance(self.file_names, list):
            files = [self.file_names]
        else:
            files = self.file_names
        # check that the number of fields match
        for file in files:
            with open(file, 'r') as f:
                line = f.readline()
                ls = line.strip().strip(self.strip_chars).split(self.delimiter)
                if self.delimiter not in line:
                    print(f"WARNING: delimiter {self.delimiter} not found in the first line of file {file}.")
                for char in self.strip_chars:
                    if char not in line:
                        print(f"WARNING: strip_char {char} not found in the first line of file {file}.")
                if len(ls) != len(self.keywords):
                    print(f"ERROR: file {file} has {len(ls)} fields, but {len(self.keywords)} keywords are provided.")
                    raise ValueError


    def read_data(self):
        """
        Function that reads the actual match data from file(s).
        """
        self.check_files()
        
        if not isinstance(self.file_names, list):
            files = [self.file_names]
        else:
            files = self.file_names

        if self.matches is not None or self.target_names is not None or self.query_names is not None:
            print("WARNING: overwriting matches.")

        self.matches = []
        self.target_names = []
        self.query_names = []
        
        for file in files:
            with open(file, 'r') as f:
                for line in f:
                    ls = line.strip().strip(self.strip_chars).split(self.delimiter)
                    match_dict = {self.keywords[i]:ls[i].strip().strip(self.strip_chars) for i in range(len(ls))}

                    # convert data to appropriate types
                    for key, val in match_dict.items():
                        match_dict[key] = standard_blast_keyword_types[key](val)

                    # create match object
                    curr_match = Match(match_dict["target"], match_dict["query"], properties=match_dict)
                    self.matches.append(curr_match)
                    self.target_names.append(match_dict["target"])
                    self.query_names.append(match_dict["query"])

        self.target_names = list(set(self.target_names))
        self.query_names = list(set(self.query_names))

        self.index_data()

    def index_data(self):
        """
        Once the data is loaded, index the data to match target names or query names
        to match objects.
        """
        if self.target_index is not None or self.query_index is not None:
            print("WARNING: overwriting indices.")
        
        self.query_index = {name:[] for name in self.query_names}
        self.target_index = {name:[] for name in self.target_names}

        for i, match in enumerate(self.matches):
            q_name = match.properties["query"]
            t_name = match.properties["target"]
            self.query_index[q_name].append(i)
            self.target_index[t_name].append(i)

    def load_sequences(self, t_file_name, q_file_name, truncate_target=True):
        """
        Given a list of target and query names, load sequences.

        t_file_name: .fasta file name for target sequences (or a list of names)
        q_file_name: .fasta file name for query sequences (or a list of names)
        """

        if not file_name_check(t_file_name) or not file_name_check(q_file_name):
            print("ERROR: problem with file names.")
            raise TypeError

        if isinstance(t_file_name, str):
            t_file_name = [t_file_name]

        if isinstance(q_file_name, str):
            q_file_name = [q_file_name]


        if self.query_names is None or self.target_names is None:
            print("ERROR: missing names. Did you read the data first?")
            raise ValueError

        if truncate_target:
            if "tstart" not in self.keywords or "tend" not in self.keywords:
                print("WARNING: no tstart/tend. Will not truncate.")
                truncate_target=False
            else:
                ranges = {}
                for t_name in self.target_names:
                    tstart = min([match.properties["tstart"] for match in self.matches if match.target == t_name])
                    tend = max([match.properties["tend"] for match in self.matches if match.target == t_name])
                    ranges[t_name] = (tstart, tend)

        print(ranges)

        self.query_sequences = {name:None for name in self.query_names}
        self.target_sequences = {name:None for name in self.target_names}

        for t_file in t_file_name:
            print(t_file)
            sset = seq_set(file_name=t_file)
            for t_name in self.target_names:
                seq = sset.find_subset_with_names([t_name], suppress_warning=True)
                if seq is None:
                    continue
                seq = sset[0]
                if truncate_target:
                    seq = seq[ranges[t_name][0]:ranges[t_name][1]]
                self.target_sequences[t_name] = seq

        for q_file in q_file_name:
            sset = seq_set(file_name=q_file)
            for q_name in self.query_names:
                seq = sset.find_subset_with_names([q_name], suppress_warning=True)
                if seq is None:
                    continue
                seq = sset[0]
                self.query_sequences[q_name] = seq

        for key, val in self.query_sequences.items():
            if val is None:
                print(f"ERROR: the sequence for query:{key} not found. Please check file names.")

        for key, val in self.target_sequences.items():
            if val is None:
                print(f"ERROR: the sequence for target:{key} not found. Please check file names.")


    def pairwise_align_matches(self, truncate_target=True):
        """
        Given matches, and their sequences, (locally) pairwise align using default blast scoring.

        If truncate_target is true, the target sequence will be truncated to tstart:tend.
        """
        if truncate_target:
            if "tstart" not in self.keywords or "tend" not in self.keywords:
                print("WARNING: tstart and/or tend not in keywoards. Target will not be truncated.")
                truncate_target = False

        for i, match in enumerate(self.matches):
            q_seq = self.query_sequences[match.query]
            t_seq = self.query_sequences[match.target]
            if truncate_target:
                tstart = match.properties['tstart']
                tend = match.properties['tend']
                t_seq = t_seq[tstart:tend]

            if q_seq.type != t_seq.type:
                print("ERROR: sequence types do not match.")
                raise TypeError

            if q_seq.type == 'dna' or q_seq.type == 'rna':
                algorithm = 'blastn'
            elif q_seq.type == 'protein':
                algorithm = 'blastp'
            else:
                print("ERROR: unknown sequence type.")
                raise TypeError

            #TODO: if frame is in keywords, then skip this?
            aligner = pairwise_alignment(t_seq, q_seq, algorithm=algorithm)
            aligner.align()
            aligner_rev = pairwise_alignment(t_seq.reverse_complement(), q_seq, algorithm=algorithm)
            aligner_rev.align()
            if aligner.score > aligner_rev.score:
                self.matches[i].alignment = copy.copy(aligner)
            else:
                self.matches[i].alignment = copy.copy(aligner)

            self.matches[i].properties[f"{algorithm}_score"] = self.matches[i].alignment.score

        #TODO: how to avoid failing the inital check, should it be run again
        self.keywords.append(f"{algorithm}_score")