"""
Genome map object that utilizes pysam and other tools.

Currently only supports SAM files but easily extendable to BAM.

TODO:
- add support for multiple ch_name/ch_nums?
- currently all the machinery assumes there is only one ch_name, but there are
no checks anywhere (which is a problem)
"""
import numpy as np
from scipy.stats import zscore
from sequence_analysis.sam_reader import SamReader
from sequence_analysis.sam_entry import SamEntry

class GenomeMap:
    """
    Class to read genome mapping info from a .sam file, check gene coverages etc.
    """
    def __init__(self, ch_name, ch_num, sample_name=None):
        self.file_names = []
        self.ch_name = ch_name
        self.ch_num = ch_num
        self.sample_name = sample_name

        self.full_lengths = None

        self.ref_len = None

        self.heatmap = None

        # range filters
        self.start = 0
        self.end = np.inf

        # minimum alignment score
        self.min_AS = 0

        # reads
        self.reads = []

        # reads that overlap with start-end are selected, but they may end up constituting a larger range
        self.read_range = (self.start, self.end)

    class Heatmap:
        """
        Inner class that makes it easier to determine heatmaps and get relevant statistics.
        """
        def __init__(self, ch_name, ch_num, start, end):
            self.ch_name = ch_name
            self.ch_num = ch_num
            self.start = start
            self.end = end

            # array to hold freq values
            # implicitly the x-values are: range(start, end)
            self.freqs = np.zeros(end - start)

        def __len__(self):
            return len(self.freqs)

        def add_read(self, read_start, read_end):
            if read_end < self.start or read_start > self.end:
                return

            self.freqs[read_start - self.start:read_end - self.start] += 1

        def get_freq_for_range(self, start, end):
            return self.freqs[start - self.start: end - self.start]

    def set_range_filter(self, start, end):
        """
        Set filter for range of nucletiodes to check.
        """
        if start > end:
            self.start = end
            self.end = start
        elif start == end:
            print("ERROR: start and end are equal to each other, nothing can be selected.")
            raise ValueError
        else:
            self.start = start
            self.end = end

        self.read_range = (self.start, self.end)

    def set_alignment_score_filter(self, min_AS):
        """
        Set filter for the minimum alignment score.

        IMPORTANT: alignment scores are assumed to be in %, but this need not be the case in any
        .sam file.
        """
        if min_AS > 100:
            print(f"WARNING: provided alignment score is {min_AS}%. Setting to 100%.")
            self.min_AS = 100
        elif min_AS < 0:
            print(f"WARNING: provided alignment score is {min_AS}%. Setting to 0%.")
            self.min_AS = 0
        else:
            self.min_AS = min_AS

    def get_reads(self, fasta_file=None, length_dict=None):
        """
        Use the sam_reader module to filter and read.
        """
        if len(self.file_names) == 0:
            print("ERROR: no input file name set.")
            raise ValueError

        self.reads = []
        headers = []

        self.read_range = (self.read_range[0], -1)

        end = self.end
        if end == np.inf:
            end = -1

        for file in self.file_names:
            sr = SamReader(file, start=self.start, end=max(end, -1), mapped_onto=self.ch_name, min_score=float(self.min_AS))
            sr.read()
            if self.min_AS > 0 and fasta_file is not None:
                sr.read_full_lengths(file_name=fasta_file)
                sr.normalize_AS_and_filter(self.min_AS)
                self.full_lengths = sr.full_lengths

            if self.min_AS > 0 and length_dict is not None:
                sr.full_lengths = length_dict
                sr.normalize_AS_and_filter(self.min_AS)
                self.full_lengths = sr.full_lengths

            self.reads += sr.sam_string_list
            self.read_range = (min(self.read_range[0], sr.seq_start), max(self.read_range[1], sr.seq_end))
            headers += sr.header

        self.reads = list(set(self.reads))
        self.end = max([int(a.split("LN:")[1]) for a in headers])
        self.ref_len = self.end

    def get_heatmap(self):
        if len(self.reads) == 0:
            print("ERROR: no reads found.")
            raise ValueError

        if self.read_range[1] == np.inf:
            raise ValueError

        self.heatmap = self.Heatmap(self.ch_name, self.ch_num, self.read_range[0], self.read_range[1])
        for read in self.reads:
            entry = SamEntry(read)
            self.heatmap.add_read(entry.pos, entry.pos + len(entry.sequence))

    def get_stats_per_range(self, ranges):
        """
        Given a list of ranges, get z-scores.
        """
        if self.heatmap is None:
            self.get_heatmap()
        averages = []
        stdevs = []
        for r in ranges:
            freq = self.heatmap.get_freq_for_range(r[0], r[1])
            avg = np.mean(freq)
            averages.append(avg)
            stdevs.append(np.std(freq))

        avg_zscore = zscore(averages)
        std_zscore = zscore(stdevs)

        return averages, avg_zscore, stdevs, std_zscore
