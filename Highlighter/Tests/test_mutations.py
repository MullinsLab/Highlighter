import unittest
from Bio import Seq, SeqRecord
from Highlighter.Mutations import Mutations, get_mutations

class MutationTests(unittest.TestCase):
    """ Test the mutation object """

    def test_mutations_start_out_empty(self):
        """ Mutations should start out empty """

        mutations = Mutations([])

        self.assertEqual(len(mutations), 0)

    def test_mutations_sequences_can_be_added(self):
        """ Mutations should be able to add sequences """

        mutations = Mutations([])
        mutations.append(SeqRecord.SeqRecord(Seq.Seq('ATGC')))

        self.assertEqual(len(mutations), 1)

    def test_get_mutations_error_if_no_sequence(self):
        """ get_mutations should error if no sequence is provided """

        with self.assertRaises(TypeError):
            get_mutations(reference='ATGC')

    def test_get_mutations_error_if_no_reference(self):
        """ get_mutations should error if no reference is provided """

        with self.assertRaises(TypeError):
            get_mutations(sequence='ATGC')

if __name__ == '__main__':
    unittest.main()