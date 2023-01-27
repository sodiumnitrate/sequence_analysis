# sequence_analysis

This python module puts together some tools to do very basic (protein/dna/rna) sequence analysis. Some of these were created for research and some for my own education. Please use at your own risk. 

## installation

    git clone https://github.com/sodiumnitrate/sequence_analysis.git
    cd sequence_analysis
    python setup.py install

## basic usage

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

## notes about performance

I have implemented global (Needleman-Wunsch) and local (Smith-Waterman) pairwise alignment algorithms in *pure python*. Their performance isn't great, as there are two nested for loops in the python code. I hope to implement these in C++ for funsies in the near future. For now, I am running Biopython's `Align.PairwiseAligner()` under the hood, using the `algorithm="biopython-global"` parameter in the `alignment` class.