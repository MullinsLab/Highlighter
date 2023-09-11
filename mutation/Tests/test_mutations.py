import unittest, os, pathlib

from itertools import product
from pprint import pprint

import mutations

from Bio import AlignIO, Phylo, SeqUtils
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Graphics import MutationPlot

def file_hash(*, file_name: str) -> str:
    """ Returns the hash of a file 
    https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file"""

    import hashlib

    hash_md5 = hashlib.md5()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    
    return hash_md5.hexdigest()

class SeqUtilsTests(unittest.TestCase):
    """ Tests of the SeqUtils module """

    def test_codon_position_error_if_invalid_sequence(self):
        """ codon_position should error if an invalid sequence is provided """

        with self.assertRaises(TypeError):
            SeqUtils.codon_position(None, 1)

    def test_codon_position_error_if_invalid_position(self):
        """ codon_position should error if an invalid position is provided """

        with self.assertRaises(TypeError):
            SeqUtils.codon_position('ATGC', None)

        with self.assertRaises(ValueError):
            SeqUtils.codon_position('ATGC', 100)

    def test_codon_position_returns_error_if_position_is_gap(self):
        """ codon_position should error if the position is a gap """

        with self.assertRaises(ValueError):
            SeqUtils.codon_position('AT-C', 2)

    def test_codon_position_returns_correct_position(self):
        """ codon_position should return the correct position """

        self.assertEqual(SeqUtils.codon_position('ATGC', 1), 1)
        self.assertEqual(SeqUtils.codon_position('A-GC', 2), 1)
        self.assertEqual(SeqUtils.codon_position('A-----GC', 7), 2)
        self.assertEqual(SeqUtils.codon_position('A-GC--A', 6),0)

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

    def test_get_mutations_returns_correct_stop_codons(self):
        """ get_mutations should return a dictionary with the correct stop codons """

        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='GTAA-', sequence='GTAA-', type='NT', stop_codons=True), {})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='GGGTAA-', sequence='GGGTAA-', type='NT', stop_codons=True), {3: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='GG--GTAA-', sequence='GG--GTAA-', type='NT', stop_codons=True), {5: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='GGGUAG-', sequence='GGGUAG-', type='NT', stop_codons=True), {3: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='TAG', sequence='TAG', type='NT', stop_codons=True), {0: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='TGA', sequence='TGA', type='NT', stop_codons=True), {0: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='TGATAG', sequence='TGATAG', type='NT', stop_codons=True), {0: ['Stop codon'], 3: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='TGATAG', sequence='TGATAG', type='NT', stop_codons=False), {})
        self.assertEqual(AlignInfo.Mutations.get_mutations(reference='TGATAG', sequence='TGATAG', type='NT'), {})

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


class MutationPlotStaticTests(unittest.TestCase):
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

        self.align_nt = AlignIO.read('mutation/Tests/Mutation/highlighter_nt.fasta', 'fasta')
        self.mutation_plot_nt = MutationPlot(self.align_nt)

        self.align_aa = AlignIO.read('mutation/Tests/Mutation/highlighter_aa.fasta', 'fasta')
        self.mutation_plot_aa = MutationPlot(self.align_aa)

        self.align_hiv_nt = AlignIO.read('mutation/Tests/Mutation/hiv_nt.fasta', 'fasta')
        self.tree_hiv_nt = Phylo.read('mutation/Tests/Mutation/hiv_nt_nexus.tre', 'nexus')
        self.mutation_plot_hiv_nt = MutationPlot(self.align_hiv_nt, tree=self.tree_hiv_nt)

        try:
            import rlPyCairo as test
            self.extended_formats: bool = True
            del test
        except ImportError:
            self.extended_formats: bool = False

    def test_mutation_plot_inits_correctly(self):
        """ MutationPlot should initialize correctly """

        self.assertIs(self.mutation_plot_nt.alignment, self.align_nt)

    def test_mutation_plot_fails_init_with_bad_tree(self):
        """ MutationPlot should fail to initialize if a bad tree is provided """

        with self.assertRaises(TypeError):
            MutationPlot(self.align_nt, tree='bad_tree')

    def test_mutation_plot_creates_valid_plot_nt(self):
        """ MutationPlot should create a valid plot """

        hashes: dict = {
            'mutation/Tests/Mutation/zz_nt.bmp': '4c560514d75c3046b55c1adefa81d91f',
            'mutation/Tests/Mutation/zz_nt.eps': '51bc03aa993d918158cd4520d4e0986a',
            'mutation/Tests/Mutation/zz_nt.gif': '081b8620d15c20ac4a9360bdd3e6b654',
            'mutation/Tests/Mutation/zz_nt.jpg': 'b8753366700dc98a2961b9d913fcb661',
            'mutation/Tests/Mutation/zz_nt.pdf': '993d78e8fcdab23345b3584f2d99c848',
            'mutation/Tests/Mutation/zz_nt.png': '9af09fdc2e6ebd3b775eda952081faab',
            'mutation/Tests/Mutation/zz_nt.ps': '51bc03aa993d918158cd4520d4e0986a',
            'mutation/Tests/Mutation/zz_nt.svg': 'a40951743cdfdd5d3b7c4d15431400fd',
            'mutation/Tests/Mutation/zz_nt.tif': 'c4856f168499d6d06d09712b6738aead',
            'mutation/Tests/Mutation/zz_nt.tiff': 'c4856f168499d6d06d09712b6738aead',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-False_stop-False.svg': 'a40951743cdfdd5d3b7c4d15431400fd',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-False_stop-True.svg': '73ce6b6b7ff466c8612bc267983cee93',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-True_stop-False.svg': 'f573fbeaa87adba49f718caf64509782',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-True_stop-True.svg': '7eaf176996901444fd6e28c6abf942a6',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-False_stop-False.svg': 'f56be954ce6a24f54dd38f3eac9d70e9',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-False_stop-True.svg': 'c5e747947afaab60d6d8b9d3206ab1cb',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-True_stop-False.svg': 'f8b0b9c2c39306ab4b07d6856094d1ea',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-True_stop-True.svg': 'b714cfc0aa70b40c19d5b6bbc533372c',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-False_stop-False.svg': 'fc498f41c78c3faf3b63d96d93fdbf5e',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-False_stop-True.svg': '39cb41c9e41e5be9497d9c0edc86f54b',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-True_stop-False.svg': '82eeb0397410605523e66a12b68747ca',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-True_stop-True.svg': 'ed5dc4da1724a5a829d1618be3e7c28f',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-False_stop-False.svg': 'f8fe33991e6641ed3af9473235901f53',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-False_stop-True.svg': 'f4124815a3f8ff5568443444621ea7ed',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-True_stop-False.svg': 'ed5f8706b458daa3231e3d66c088c1c1',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-True_stop-True.svg': '05ec74319594752382635eb819c61e1a'
        }

        # PDFs have a creation date in them, so can't be checked by hash as they're different each time
        formats = ["PS", "EPS", "SVG", "JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"] # , "PDF"

        for format in formats:
            if not self.extended_formats and format in ["JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"]:
                continue

            with self.subTest(format=format):
                file_name: str = f"mutation/Tests/Mutation/zz_nt.{format.lower()}"
                self.mutation_plot_nt.draw(file_name, output_format=format)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        titles: list = [None, "A Title"]
        apobecs: list = [True, False]
        g_to_as: list = [True, False]
        stop_codonss: list = [True, False]

        for title, apobec, g_to_a, stop_codons in product(titles, apobecs, g_to_as, stop_codonss):
            with self.subTest(title=title, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons):
                file_name: str = f"mutation/Tests/Mutation/zz_nt_title-{bool(title)}_apobec-{apobec}_gtoa-{g_to_a}_stop-{stop_codons}.svg"
                self.mutation_plot_nt.draw(file_name, title=title, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

    def test_mutation_plot_creates_valid_plot_aa(self):
        """ MutationPlot should create a valid plot """

        hashes: dict = {
            'mutation/Tests/Mutation/zz_aa.bmp': '6202b8038419e4845f56c6f50dfc6d12',
            'mutation/Tests/Mutation/zz_aa.eps': 'a626c585761825a7c42480fb0f89d7cc',
            'mutation/Tests/Mutation/zz_aa.gif': '8890a7db5610b0afc3120dd71de8c9e0',
            'mutation/Tests/Mutation/zz_aa.jpg': '63163e9f52a96a61eaaffaa1323e5bff',
            'mutation/Tests/Mutation/zz_aa.png': '14e5519b9d1ed682edadb5d915950b7a',
            'mutation/Tests/Mutation/zz_aa.ps': 'a626c585761825a7c42480fb0f89d7cc',
            'mutation/Tests/Mutation/zz_aa.svg': 'c3ecb8cd9b76f70631d9ac0e78318ac3',
            'mutation/Tests/Mutation/zz_aa.tif': '5c94dad4c1608207aa627b1af40078ab',
            'mutation/Tests/Mutation/zz_aa.tiff': '5c94dad4c1608207aa627b1af40078ab',
            'mutation/Tests/Mutation/zz_aa_title-False_glycosylation-False.svg': 'c3ecb8cd9b76f70631d9ac0e78318ac3',
            'mutation/Tests/Mutation/zz_aa_title-False_glycosylation-True.svg': '6baa39bb672beef0a579180a21fddc7f',
            'mutation/Tests/Mutation/zz_aa_title-True_glycosylation-False.svg': '223bbaff099f3267481c1b5076f0cf63',
            'mutation/Tests/Mutation/zz_aa_title-True_glycosylation-True.svg': 'ee6576d9f7ef8cdc908507ae22df4d8a'
        }

        # PDFs have a creation date in them, so can't be checked by hash as they're different each time
        formats = ["PS", "EPS", "SVG", "JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"] # , "PDF"

        for format in formats:
            if not self.extended_formats and format in ["JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"]:
                continue

            with self.subTest(format=format):
                file_name: str = f"mutation/Tests/Mutation/zz_aa.{format.lower()}"
                self.mutation_plot_aa.draw(file_name, output_format=format)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        titles: list = [None, "A Title"]
        glycosylations: list = [True, False]

        for title, glycosylation in product(titles, glycosylations):
            with self.subTest(title=title, glycosylation=glycosylation):
                file_name: str = f"mutation/Tests/Mutation/zz_aa_title-{bool(title)}_glycosylation-{glycosylation}.svg"
                self.mutation_plot_aa.draw(file_name, title=title, glycosylation=glycosylation)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

if __name__ == '__main__':
    unittest.main()