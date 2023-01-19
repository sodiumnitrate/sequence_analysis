# sequence_analysis

This python module puts together some tools to do very basic (protein/dna/rna) sequence analysis. Some of these were created for research and some for my own education. Please use at your own risk. :upside_down_face:

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