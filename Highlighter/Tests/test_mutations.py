import unittest
import pathlib

import Mutations

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from Highlighter.Mutations import get_mutations

class MutationStaticTests(unittest.TestCase):
    """ Test the static methods of the mutation object """

    def test_get_mutations_error_if_no_sequence(self):
        """ get_mutations should error if no sequence is provided """

        with self.assertRaises(TypeError):
            Mutations.get_mutations(reference='ATGC')

    def test_get_mutations_error_if_no_reference(self):
        """ get_mutations should error if no reference is provided """

        with self.assertRaises(TypeError):
            Mutations.get_mutations(sequence='ATGC')

    def test_get_mutations_error_if_bad_types(self):
        """ get_mutations should error if bad types are provided """

        with self.assertRaises(TypeError):
            Mutations.get_mutations(sequence=1, reference='ATGC')

        with self.assertRaises(TypeError):
            Mutations.get_mutations(sequence='ATGC', reference=[])

    def test_get_mutations_error_if_different_lengths(self):
        """ get_mutations should error if the sequence and reference are different lengths """

        with self.assertRaises(ValueError):
            Mutations.get_mutations(sequence='ATGC', reference='ATG')

    def test_get_mutations_succeeds_with_mixed_types(self):
        """ get_mutations should succeed with mixed types """

        try:
            Mutations.get_mutations(sequence='ATGC', reference=Seq('ATGC'))
        except:
            self.fail('get_mutations failed with mixed types')

        try:
            Mutations.get_mutations(sequence=SeqRecord(Seq('ATGC')), reference='ATGC')
        except:
            self.fail('get_mutations failed with mixed types')

    def test_get_mutations_returns_dict(self):
        """ get_mutations should return a dictionary """

        self.assertEqual(type(Mutations.get_mutations(sequence='ATGC', reference='ATGC')), dict)

    def test_get_mutations_returns_empty_dict_given_same_sequence_and_reference(self):
        """ get_mutations should return an empty dictionary if the sequence and reference are the same """

        self.assertEqual(Mutations.get_mutations(sequence='ATGC', reference='ATGC'), {})

    def test_get_mutations_returns_non_empty_dict_given_different_sequence_and_reference(self):
        """ get_mutations should return a non-empty dictionary if the sequence and reference are different """

        self.assertNotEqual(Mutations.get_mutations(sequence='ATGC', reference='ATGG'), {})

    def test_get_mutations_returns_correct_dict_given_different_sequence_and_reference(self):
        """ get_mutations should return a dictionary with the correct keys and values """

        self.assertEqual(Mutations.get_mutations(reference='GTGCGGC-', sequence='AATGCA-T'), {0: ['A'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(Mutations.get_mutations(g_to_a=True, reference='GTGCGGC-', sequence='AATGCA-T'), {0: ['A', 'G->A mutation'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(Mutations.get_mutations(apobec=True, reference='GTGCGGC-', sequence='AATGCA-T'), {0: ['A', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(Mutations.get_mutations(g_to_a=True, apobec=True, reference='GTGCGGC-', sequence='AATGCA-T'), {0: ['A', 'G->A mutation', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})

    # def test_stuff_with_bio_align(self):
    #     """ Test stuff with Bio.Align """

    #     align = AlignIO.read(pathlib.PurePath(pathlib.Path(__file__).parent.resolve(), 'Mutations/test.fasta'), 'fasta')

    #     self.assertEqual(align.prep_mutations(), "test")

if __name__ == '__main__':
    unittest.main()