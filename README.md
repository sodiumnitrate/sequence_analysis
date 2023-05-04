# sequence_analysis

This python module puts together some tools to do very basic (protein/dna/rna) sequence analysis. Some of these were created for research and some for my own education. Please use at your own risk. 

## Installation

    git clone https://github.com/sodiumnitrate/sequence_analysis.git
    cd sequence_analysis
    python setup.py install

## Basic usage

Loading module

    import sequence_analysis as sa

Reading in a `.fasta` file named `test.fasta`

    sequences = sa.seq_set.seq_set(file_name="test.fasta")

The type of the sequence should be automatically detected, and can be checked by:

    sequences.type

Creating a sequence object from a string (e.g. `"ALKALI"`)

    seq = sa.sequence.sequence("ALKALI")

Writing a `seq_set` object to `output.fasta`

    sequences.write_fasta("output.fasta")

### Finding open reading frames

I have created an `OpenReadingFrame` class find and store a valid open reading frame, storing information such as: `strand`, `frame`, as well as the starting and ending indices of the codons for a given frame, and the corresponding protein sequence. The `sequence` class contains a method, `detailed_orf_finder()` to find all open reading frames. The following example illustrates how to use it. We will fetch the nucleotide sequence of a human ubiquitin from GenBank, and find its open reading frames.

First, import the necessary modules

    from sequence_analysis.sequence import sequence
    from sequence_analysis.open_reading_frame import OpenReadingFrame
    from sequence_analysis.genbank_entry import GenBankEntry

Then, fetch the GenBank info

    ubb = GenBankEntry('NM_001281716.2')
    ubb.fetch("your.email@provider.com", skip_origin=False)
    dna_seq = ubb.dna_sequence

Find open reading frames. Note that we are setting `min_orf_len=80` to reproduce the results of NCBI's ORFfinder:

    orfs = dna_seq.detailed_orf_finder(min_orf_len=80)

This results in 6 open reading frames:

    >>> len(orfs)
    6

This is consistent with the results of NCBI's server. Print results by doing

    >>> for i, orf in enumerate(orfs):
    ...     print(i, orf.protein_sequence)
    ...
    0 MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGC*
    1 MSLCWSGGMPSLSWILAFTFSMVSLGSTSRVMVLPVRVFTKICIPPLRRRTRCRVDSFWML*
    2 MSLCWSGGMPSLSWILAFTFSMVSLGSTSRVMVLPVRVFTKICIPPLRRRTRCRVDSFWML*
    3 MSLCWSGGMPSLSWILAFTFSMVSLGSTSRVMVLPVRVFTKICILTPHVDAVSVRR*
    4 MAFAVPSDGITLHYSHLPQLKFRNYKFQ*
    5 MPSLGTANAMTEELTATPQAQDQVQGRLFLDVVVRKSAAIFQLLACKDEPLLVGGDAFFILDLGLHIFDGVTGLHLQSDGLAGQGLHEDLHTTSQTQDQVQGRLLLDVVVRKSAAIFQLLACKDEPLLVGRDAFFILDLGLHIFDGVTGLHFQGDGLAGQGLHEDLHTTSQTQDQVQGRLLLDVVVRKSTAIFQLLACKDEPLLVGGNAFLILDLGLHIFDGVTGLHLKGDGLAGKGFHEDLHFDPSRRRRLRAPLKVVRRAQPPERQFRLFN*

These are indeed what NCBI ORFfinder finds. You can also check and compare the strand, frame, start, and stop, for the first ORF we find, by

    orfs[0].strand
    orfs[0].frame
    orfs[0].start
    orfs[0].stop

## Notes about performance

I have implemented global (Needleman-Wunsch) and local (Smith-Waterman) pairwise alignment algorithms in *pure python*. Their performance isn't great, as there are two nested for loops in the python code. I hope to implement these in C++ for funsies in the near future. For now, I am running Biopython's `Align.PairwiseAligner()` under the hood, using the `algorithm="biopython-global"` parameter in the `alignment` class.