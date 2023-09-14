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

class MutationStaticMismatchTests(unittest.TestCase):
    """ Test the static mismatch methods of the mutation object """

    def test_get_mismatches_error_if_no_type_given(self):
        """ get_mismatches should error if no type is provided """

        with self.assertRaises(ValueError):
            AlignInfo.Mutations.get_mismatches(sequence='ATGC', references=['ATGC'])

    def test_get_mismatches_error_if_no_sequence(self):
        """ get_mismatches should error if no sequence is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mismatches(references='ATGC', seq_type='NT')

    def test_get_mismatches_error_if_no_reference(self):
        """ get_mismatches should error if no reference is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mismatches(sequence='ATGC', seq_type='NT')

    def test_get_mismatches_error_if_bad_types(self):
        """ get_mismatches should error if bad types are provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mismatches(sequence=1, references='ATGC', seq_type='NT')

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_mismatches(sequence='ATGC', references=[], seq_type='NT')

    def test_get_mismatches_error_if_different_lengths(self):
        """ get_mismatches should error if the sequence and reference are different lengths """

        with self.assertRaises(ValueError):
            AlignInfo.Mutations.get_mismatches(sequence='ATGC', references='ATG', seq_type='NT')

    def test_get_mismatches_succeeds_with_mixed_types(self):
        """ get_mismatches should succeed with mixed types """

        try:
            AlignInfo.Mutations.get_mismatches(sequence='ATGC', references=Seq('ATGC'), seq_type='NT')
        except:
            self.fail('get_mismatches failed with mixed types')

        try:
            AlignInfo.Mutations.get_mismatches(sequence=SeqRecord(Seq('ATGC')), references='ATGC', seq_type='NT')
        except:
            self.fail('get_mismatches failed with mixed types')

    def test_get_mismatches_returns_dict(self):
        """ get_mismatches should return a dictionary """

        self.assertEqual(type(AlignInfo.Mutations.get_mismatches(sequence='ATGC', references='ATGC', seq_type='NT')), dict)

    def test_get_mismatches_returns_empty_dict_given_same_sequence_and_reference(self):
        """ get_mismatches should return an empty dictionary if the sequence and reference are the same """

        self.assertEqual(AlignInfo.Mutations.get_mismatches(sequence='ATGC', references='ATGC', seq_type='NT'), {})

    def test_get_mismatches_returns_non_empty_dict_given_different_sequence_and_reference(self):
        """ get_mismatches should return a non-empty dictionary if the sequence and reference are different """

        self.assertNotEqual(AlignInfo.Mutations.get_mismatches(sequence='ATGC', references='ATGG', seq_type='NT'), {})

    def test_get_mismatches_returns_correct_glycosylation_sites(self):
        """ get_mismatches should return a dictionary with the correct glycosylation sites """

        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='GNS-SQ', sequence='GNS-SQ', seq_type='AA', glycosylation=True), {1: ['Glycosylation']})

    def test_get_mismatches_returns_correct_stop_codons(self):
        """ get_mismatches should return a dictionary with the correct stop codons """

        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='GTAA-', sequence='GTAA-', seq_type='NT', stop_codons=True), {})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='GGGTAA-', sequence='GGGTAA-', seq_type='NT', stop_codons=True), {3: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='GG--GTAA-', sequence='GG--GTAA-', seq_type='NT', stop_codons=True), {5: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='GGGUAG-', sequence='GGGUAG-', seq_type='NT', stop_codons=True), {3: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='TAG', sequence='TAG', seq_type='NT', stop_codons=True), {0: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='TGA', sequence='TGA', seq_type='NT', stop_codons=True), {0: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='TGATAG', sequence='TGATAG', seq_type='NT', stop_codons=True), {0: ['Stop codon'], 3: ['Stop codon']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='TGATAG', sequence='TGATAG', seq_type='NT', stop_codons=False), {})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='TGATAG', sequence='TGATAG', seq_type='NT'), {})

    def test_get_mismatches_returns_correct_dict_given_different_sequence_and_reference_nt(self):
        """ get_mismatches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(g_to_a=True, references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A', 'G->A mutation'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(apobec=True, references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Mutations.get_mismatches(g_to_a=True, apobec=True, references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A', 'G->A mutation', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})

    def test_get_mismatches_returns_correct_dict_given_different_sequence_and_reference_aa(self):
        """ get_mismatches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_mismatches(references='MRVMEIRRNYQHL--', sequence='MRAMK-RRNYQHL--', seq_type='AA'), {2: ['A'], 4: ['K'], 5: ['Gap']})


class MutationStaticMatchTests(unittest.TestCase):
    """ Test the static match methods of the mutation object """

    def test_get_matches_error_if_no_type_given(self):
        """ get_matches should error if no type is provided """

        with self.assertRaises(ValueError):
            AlignInfo.Mutations.get_matches(sequence='ATGC', references='ATGC')

    def test_get_matches_error_if_no_sequence(self):
        """ get_matches should error if no sequence is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_matches(references='ATGC', seq_type='NT')

    def test_get_matches_error_if_no_reference(self):
        """ get_matches should error if no reference is provided """
    
        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_matches(sequence='ATGC', seq_type='NT')

    def test_get_matches_error_if_bad_types(self):
        """ get_matches should error if bad types are provided """

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_matches(sequence=1, references='ATGC', seq_type='NT')

        with self.assertRaises(TypeError):
            AlignInfo.Mutations.get_matches(sequence='ATGC', references=1, seq_type='NT')

    def test_get_matches_error_if_different_lengths(self):
        """ get_matches should error if the sequence and reference are different lengths """

        with self.assertRaises(ValueError):
            AlignInfo.Mutations.get_matches(sequence='ATGC', references='ATG', seq_type='NT')

    def test_get_matches_succeeds_with_mixed_types(self):
        """ get_matches should succeed with mixed types """

        try:
            AlignInfo.Mutations.get_matches(sequence='ATGC', references=Seq('ATGC'), seq_type='NT')
        except:
            self.fail('get_matches failed with mixed types')

        try:
            AlignInfo.Mutations.get_matches(sequence=SeqRecord(Seq('ATGC')), references='ATGC', seq_type='NT')
        except:
            self.fail('get_matches failed with mixed types')

    def test_get_matches_returns_dict(self):
        """ get_matches should return a dictionary """

        self.assertEqual(type(AlignInfo.Mutations.get_matches(sequence='ATGC', references='ATGC', seq_type='NT')), dict)

    def test_get_matches_returns_non_empty_dict_given_different_sequence_and_reference(self):
        """ get_mismatches should return a non-empty dictionary if the sequence and reference are different """

        self.assertNotEqual(AlignInfo.Mutations.get_matches(sequence='ATGC', references='ATGG', seq_type='NT'), {})

    def test_get_matches_returns_correct_matches_given_different_sequence_and_reference(self):
        """ get_matches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_matches(references='GTGCGGC-', sequence='GTGCGGCT', seq_type='NT'), {7: [0]})

    def test_get_matches_returns_correct_matches_given_multiple_references(self):
        """ get_matches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_matches(references=['GTGCGGC-', 'GTGTGGCT'], sequence='GTGCCGCT', seq_type='NT'), {3: [1], 4: [0, 1], 7: [0]})

   
class MutationObjectMismatchTests(unittest.TestCase):
    """ Tests that use the mutation object """

    def setUp(self):
        """ Set up an align object to use for testing """

        self.short_align = AlignIO.read(pathlib.PurePath(pathlib.Path(__file__).parent.resolve(), 'Mutation/short_test_nt.fasta'), 'fasta')
        self.short_mutations = AlignInfo.Mutations(self.short_align, seq_type='NT')

        self.medium_align = AlignIO.read(pathlib.PurePath(pathlib.Path(__file__).parent.resolve(), 'Mutation/medium_test_nt.fasta'), 'fasta')
        self.medium_mutations = AlignInfo.Mutations(self.medium_align, seq_type='NT')

    def test_list_mismatches_raises_error_on_bad_reference(self):
        """ list_mismatches should raise an error if the reference is not in the alignment """

        with self.assertRaises(IndexError):
            self.short_mutations.list_mismatches(references='bad_test')

        with self.assertRaises(IndexError):
            self.short_mutations.list_mismatches(references=7)

    def test_list_matches_raises_error_on_bad_reference(self):
        """ list_matches should raise an error if the reference is not in the alignment """

        with self.assertRaises(IndexError):
            self.short_mutations.list_matches(references='bad_test')

        with self.assertRaises(IndexError):
            self.short_mutations.list_matches(references=7)

    def test_list_mismatches_returns_correct_list(self):
        """ list_mismatches should return a the correct list """

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

        for reference, apobec, g_to_a in product([0, 1], [False, True], [False, True]):
            with self.subTest(references=reference, apobec=apobec, g_to_a=g_to_a):
                self.assertEqual(self.short_mutations.list_mismatches(references=reference, apobec=apobec, g_to_a=g_to_a), result[reference][apobec][g_to_a])

    def test_list_matches_returns_correct_list(self):
        """ list_matches should return a the correct list """

        result: list = [{}, {}, {0: [1], 7: [0]}, {0: [0], 7: [1]}, {0: [0], 3: [0, 1], 7: [1]}]
        self.assertEqual(self.medium_mutations.list_matches(references=[0, 1]), result)

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
            from reportlab.graphics import renderPM as test
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

        #return

        hashes: dict = {
            'mutation/Tests/Mutation/zz_nt.bmp': '81c8c577e086484c01094058eabfca59',
            'mutation/Tests/Mutation/zz_nt.eps': '51bc03aa993d918158cd4520d4e0986a',
            'mutation/Tests/Mutation/zz_nt.gif': 'def41a67d6d8314892cd4b8c61fdf7c1',
            'mutation/Tests/Mutation/zz_nt.jpg': '37caaa5fd4abdfe10ff1c3c45176626f',
            'mutation/Tests/Mutation/zz_nt.png': 'ca63bd39385205c7b03b120e5fde50ee',
            'mutation/Tests/Mutation/zz_nt.ps': '51bc03aa993d918158cd4520d4e0986a',
            'mutation/Tests/Mutation/zz_nt.svg': 'a40951743cdfdd5d3b7c4d15431400fd',
            'mutation/Tests/Mutation/zz_nt.tif': 'b42ed4b0f0e7b28fedd57d2f4a5d3f91',
            'mutation/Tests/Mutation/zz_nt.tiff': '105b7c4e7822f5b40df1510fd6dfe199',
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
                self.mutation_plot_nt.draw_mismatches(file_name, output_format=format)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        titles: list = [None, "A Title"]
        apobecs: list = [True, False]
        g_to_as: list = [True, False]
        stop_codonss: list = [True, False]

        for title, apobec, g_to_a, stop_codons in product(titles, apobecs, g_to_as, stop_codonss):
            with self.subTest(title=title, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons):
                file_name: str = f"mutation/Tests/Mutation/zz_nt_title-{bool(title)}_apobec-{apobec}_gtoa-{g_to_a}_stop-{stop_codons}.svg"
                self.mutation_plot_nt.draw_mismatches(file_name, title=title, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)


    def test_mutation_plot_creates_valid_plot_aa(self):
        """ MutationPlot should create a valid plot """

        #return

        hashes: dict = {
            'mutation/Tests/Mutation/zz_aa.bmp': '1a203b5b4ff33d19dea7565230d45a03',
            'mutation/Tests/Mutation/zz_aa.eps': 'a626c585761825a7c42480fb0f89d7cc',
            'mutation/Tests/Mutation/zz_aa.gif': '91db6b2334176c1194e16c11f9e5d950',
            'mutation/Tests/Mutation/zz_aa.jpg': 'ee9c52e982400d03efa1e27b9fab3294',
            'mutation/Tests/Mutation/zz_aa.png': 'afa2d62a43bba2c75fac67702752c5f7',
            'mutation/Tests/Mutation/zz_aa.ps': 'a626c585761825a7c42480fb0f89d7cc',
            'mutation/Tests/Mutation/zz_aa.svg': 'c3ecb8cd9b76f70631d9ac0e78318ac3',
            'mutation/Tests/Mutation/zz_aa.tif': '15daa35eaa767dc1f2d055a0d04679c8',
            'mutation/Tests/Mutation/zz_aa.tiff': '63103557e2ae5acf9959dc00e3dcb440',
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
                self.mutation_plot_aa.draw_mismatches(file_name, output_format=format)
                hashes[file_name] = file_hash(file_name=file_name)
                os.remove(file_name)

        titles: list = [None, "A Title"]
        glycosylations: list = [True, False]

        for title, glycosylation in product(titles, glycosylations):
            with self.subTest(title=title, glycosylation=glycosylation):
                file_name: str = f"mutation/Tests/Mutation/zz_aa_title-{bool(title)}_glycosylation-{glycosylation}.svg"
                self.mutation_plot_aa.draw_mismatches(file_name, title=title, glycosylation=glycosylation)
                hashes[file_name] = file_hash(file_name=file_name)
                os.remove(file_name)


if __name__ == '__main__':
    unittest.main()