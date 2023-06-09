"""
This file contains unit tests for the GibbsSampler class and its methods.
"""

from sequence_analysis.gibbs_sampler import GibbsSampler
from sequence_analysis.seq_set import seq_set
from sequence_analysis.sequence import sequence

import pdb

class TestGibbsSampler:
    def test_init(self):
        sequences = [sequence('AAAAAAAAALKALIAAAAA'), sequence('MMMMMMMMMALKALIMMMMM'),
                     sequence('EEALKALIEEEEEEEEE'), sequence('WALKALIWWIIWWW')]
        sequences = seq_set(sequences)

        sampler = GibbsSampler(sequences, 6)

        assert sequences.type == 'protein'
        assert len(sampler.alphabet) == 20

    def test_initialize_sampler(self):
        sequences = [sequence('AAAAAAAAALKALIAAAAA'), sequence('MMMMMMMMMALKALIMMMMM'),
                     sequence('EEALKALIEEEEEEEEE'), sequence('WALKALIWWIIWWW')]
        sequences = seq_set(sequences)

        sampler = GibbsSampler(sequences, 6)
        sampler.initialize_sampler(random_init=False)
        assert sampler.indices is None

        sampler.initialize_sampler()

        assert len(sampler.indices) == len(sequences)

        sampler.initialize_sampler(random_init=False, indices=[1,1,1,1])
        assert sampler.indices[0] == 1
        assert sampler.indices[3] == 1

    def test_calculate_background_probabilities(self):
        sequences = [sequence('AAAAAAAAALKALIAAAAA'), sequence('MMMMMMMMMALKALIMMMMM'),
                     sequence('EEALKALIEEEEEEEEE'), sequence('WALKALIWWIIWWW')]
        sequences = seq_set(sequences)

        sampler = GibbsSampler(sequences, 6)
        sampler.initialize_sampler()

        sampler.calculate_background_probabilities(list(range(len(sequences))))

        assert abs(sum(sampler.background_probabilities)-1) < 1e-10

    def test_calculate_pssm(self):
        sequences = [sequence('AAAAAAAAALKALIAAAAA'), sequence('MMMMMMMMMALKALIMMMMM'),
                     sequence('EEALKALIEEEEEEEEE'), sequence('WALKALIWWIIWWW')]
        sequences = seq_set(sequences)

        sampler = GibbsSampler(sequences, 6)
        sampler.initialize_sampler()
        sampler.calculate_pssm()

        assert len(sampler.pssm) == len(sampler.alphabet)
        assert len(sampler.pssm[0]) == 6

        # check values
        column_sums = [sum(l[i]* sampler.background_probabilities[j] for j,l in enumerate(sampler.pssm)) for i in range(6)]
        for val in column_sums:
            assert abs(val - 1) < 1e-10

    def test_calculate_score(self):
        sequences = [sequence('AAAAAAAAALKALIAAAAA'), sequence('MMMMMMMMMALKALIMMMMM'),
                     sequence('EEALKALIEEEEEEEEE'), sequence('WALKALIWWIIWWW')]
        sequences = seq_set(sequences)

        sampler = GibbsSampler(sequences, 6)
        sampler.initialize_sampler()
        sampler.calculate_pssm()

        sampler.calculate_overall_score()

        assert sampler.total_score is not None

    def test_iterate(self):
        sequences = [sequence('AAAAAAAAALKALIAAAAA'), sequence('MMMMMMMMMALKALIMMMMM'),
                     sequence('EEALKALIEEEEEEEEE'), sequence('WALKALIWWIIWWW')]
        sequences = seq_set(sequences)

        sampler = GibbsSampler(sequences, 6)
        sampler.initialize_sampler()
        sampler.calculate_pssm()

        sampler.iterate()

        assert sampler.converged