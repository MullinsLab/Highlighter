import unittest
import pathlib

import mutations

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Graphics import MutationPlot

class MutationStaticTests(unittest.TestCase):
    """ Test the static methods of the mutation object """

    def test_get_mutations_error_if_no_type_given(self):
        """ get_mutations should error if no type is provided """

        with self.assertRaises(ValueError):
            AlignInfo.Mutations.get_mutations(sequence='ATGC', reference='ATGC')

    def test_get_mutations_error_if_no_sequence(self):
        """ get_mutations should error if no sequence is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mutations(reference='ATGC', type='NT')

    def test_get_mutations_error_if_no_reference(self):
        """ get_mutations should error if no reference is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mutations(sequence='ATGC', type='NT')

    def test_get_mutations_error_if_bad_types(self):
        """ get_mutations should error if bad types are provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mutations(sequence=1, reference='ATGC', type='NT')

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mutations(sequence='ATGC', reference=[], type='NT')

    def test_get_mutations_error_if_different_lengths(self):
        """ get_mutations should error if the sequence and reference are different lengths """

        with self.assertRaises(ValueError):
            AlignInfo.Mutations.get_mutations(sequence='ATGC', reference='ATG', type='NT')

    def test_get_mutations_succeeds_with_mixed_types(self):
        """ get_mutations should succeed with mixed types """

        try:
            AlignInfo.Mutations.get_mutations(sequence='ATGC', reference=Seq('ATGC'), type='NT')
        except:
            self.fail('get_mutations failed with mixed types')

        try:
            AlignInfo.Mutations.get_mutations(sequence=SeqRecord(Seq('ATGC')), reference='ATGC', type='NT')
        except:
            self.fail('get_mutations failed with mixed types')

    def test_get_mutations_returns_dict(self):
        """ get_mutations should return a dictionary """

        self.assertEqual(type(AlignInfo.Mutations.get_mutations(sequence='ATGC', reference='ATGC', type='NT')), dict)

    def test_get_mutations_returns_empty_dict_given_same_sequence_and_reference(self):
        """ get_mutations should return an empty dictionary if the sequence and reference are the same """

        self.assertEqual(AlignInfo.Mutations.get_mutations(sequence='ATGC', reference='ATGC', type='NT'), {})

    def test_get_mutations_returns_non_empty_dict_given_different_sequence_and_reference(self):
        """ get_mutations should return a non-empty dictionary if the sequence and reference are different """

        self.assertNotEqual(AlignInfo.Mutations.get_mutations(sequence='ATGC', reference='ATGG', type='NT'), {})

    def test_get_mutations_returns_correct_glycosylation_sites(self):
        """ get_mutations should return a dictionary with the correct glycosylation sites """

        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='GNS-SQ', sequence='GNS-SQ', type='AA', glycosylation=True), {1: ['Glycosylation']})

    def test_get_mutations_returns_correct_dict_given_different_sequence_and_reference_nt(self):
        """ get_mutations should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='GTGCGGC-', sequence='AATGCA-T', type='NT'), {0: ['A'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(g_to_a=True, reference='GTGCGGC-', sequence='AATGCA-T', type='NT'), {0: ['A', 'G->A mutation'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(apobec=True, reference='GTGCGGC-', sequence='AATGCA-T', type='NT'), {0: ['A', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(g_to_a=True, apobec=True, reference='GTGCGGC-', sequence='AATGCA-T', type='NT'), {0: ['A', 'G->A mutation', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})

    def test_get_mutations_returns_correct_dict_given_different_sequence_and_reference_aa(self):
        """ get_mutations should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='MRVMEIRRNYQHL--', sequence='MRAMK-RRNYQHL--', type='AA'), {2: ['A'], 4: ['K'], 5: ['Gap']})
        
class MutationObjectTests(unittest.TestCase):
    """ Tests that use the mutation object """

    def setUp(self):
        """ Set up an align object to use for testing """

        self.short_align = AlignIO.read(pathlib.PurePath(pathlib.Path(__file__).parent.resolve(), 'Mutation/short_test_nt.fasta'), 'fasta')
        self.short_mutations = AlignInfo.Mutations(self.short_align, type='NT')

    def test_list_mutations_raises_error_on_bad_reference(self):
        """ list_mutations should raise an error if the reference is not in the alignment """

        with self.assertRaises(IndexError):
            self.short_mutations.list_mutations(reference='bad_test')

        with self.assertRaises(IndexError):
            self.short_mutations.list_mutations(reference=7)

    def test_list_mutations_returns_correct_list(self):
        """ list_mutations should return a list of the correct length """

        result: dict = {
            0: {
                False: {
                    False: [{}, {0: ['A'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T'], 11: ["C"]}],
                    True: [{}, {0: ['A', 'G->A mutation'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T'], 11: ["C"]}]
                },
                True: {
                    False: [{}, {0: ['A', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T'], 11: ["C"]}],
                    True: [{}, {0: ['A', 'G->A mutation', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T'], 11: ["C"]}]
                },
            },
            1: {
                False: {
                    False: [{0: ['G'], 1: ['T'], 2: ['G'], 3: ['C'], 4: ['G'], 5: ['G'], 6: ['C'], 7: ['Gap'], 11: ["T"]}, {}],
                    True: [{0: ['G'], 1: ['T'], 2: ['G'], 3: ['C'], 4: ['G'], 5: ['G'], 6: ['C'], 7: ['Gap'], 11: ["T"]}, {}]
                },
                True: {
                    False: [{0: ['G'], 1: ['T'], 2: ['G'], 3: ['C'], 4: ['G'], 5: ['G'], 6: ['C'], 7: ['Gap'], 11: ["T"]}, {}],
                    True: [{0: ['G'], 1: ['T'], 2: ['G'], 3: ['C'], 4: ['G'], 5: ['G'], 6: ['C'], 7: ['Gap'], 11: ["T"]}, {}]
                }
            }
        }

        for reference in [0, 1]:
            for apobec in [False, True]:
                for g_to_a in [False, True]:
                    with self.subTest(reference=reference, apobec=apobec, g_to_a=g_to_a):
                        self.assertEqual(self.short_mutations.list_mutations(reference=reference, apobec=apobec, g_to_a=g_to_a), result[reference][apobec][g_to_a])


class MutationPlotStateicTests(unittest.TestCase):
    """ Tests that use the static methods of the MutationPlot class """

    def test_mutation_plot_guesses_alignment_type_correctly(self):
        """ MutationPlot guesses the alignment type correctly """

        self.assertEqual(MutationPlot.guess_alignment_type(AlignIO.read('mutation/Tests/Mutation/highlighter_nt.fasta', 'fasta')), 'NT')
        self.assertEqual(MutationPlot.guess_alignment_type(AlignIO.read('mutation/Tests/Mutation/highlighter_aa.fasta', 'fasta')), 'AA')

    def test_mutation_plot_significant_digits(self):
        """ MutationPlot significant_digits returns correctly """

        self.assertEqual(MutationPlot.significant_digits(1234), 1200)
        self.assertEqual(MutationPlot.significant_digits(51234), 51000)


class MutationPlotTests(unittest.TestCase):
    """ Holds test of the MutationPlot class """

    def setUp(self):
        """ Set up an align object to use for testing """

        self.align = AlignIO.read('mutation/Tests/Mutation/test.fasta', 'fasta')
        self.mutation_plot = MutationPlot(self.align)

    def test_mutation_plot_inits_correctly(self):
        """ MutationPlot should initialize correctly """

        self.assertIs(self.mutation_plot.alignment, self.align)

    def test_mutation_plot_fails_init_with_bad_tree(self):
        """ MutationPlot should fail to initialize if a bad tree is provided """

        with self.assertRaises(TypeError):
            MutationPlot(self.align, tree='bad_tree')

    # def test_mutation_plot_height_is_18_per_sequence(self):
    #     """ MutationPlot should have a height of 18 per sequence """

    #     self.assertEqual(self.mutation_plot.height, 18 * len(self.align))

if __name__ == '__main__':
    unittest.main()