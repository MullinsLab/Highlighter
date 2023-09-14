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

        self.assertEqual(AlignInfo.Mutations.get_matches(references='GTGCGGC-', sequence='GTGCGGCT', seq_type='NT'), {7: [0, -1]})

    def test_get_matches_returns_correct_matches_given_multiple_references(self):
        """ get_matches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Mutations.get_matches(references=['GTGCGGC-', 'GTGTGGCT'], sequence='GTGCCGCT', seq_type='NT'), {3: [1], 4: [0, 1, -1], 7: [0]})

   
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

        result: list = [{}, {}, {0: [1], 7: [0]}, {0: [0], 7: [1]}, {0: [0], 3: [0, 1, -1], 7: [1]}]
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

        # return

        hashes: dict = {
            'mutation/Tests/Mutation/zz_nt.bmp': '151197a072cf41782f7a563b9c7835fa',
            'mutation/Tests/Mutation/zz_nt.eps': '5673fef3bdd0121d0a5a5e9d84746e7b',
            'mutation/Tests/Mutation/zz_nt.gif': 'ffd504b515e8802be0421accdab17663',
            'mutation/Tests/Mutation/zz_nt.jpg': '835bb0df30ac33166516835e93d8093d',
            'mutation/Tests/Mutation/zz_nt.png': 'e9ad40ff6000198b1aa548dd89c8cc4e',
            'mutation/Tests/Mutation/zz_nt.ps': '5673fef3bdd0121d0a5a5e9d84746e7b',
            'mutation/Tests/Mutation/zz_nt.svg': '89399497a7a9ccf4a3bed6c5e9f62a28',
            'mutation/Tests/Mutation/zz_nt.tif': 'd804ca950e045af0c84ff852bcd2c9bb',
            'mutation/Tests/Mutation/zz_nt.tiff': '1c9aeb61eb111ffd7692a05d22b7bbb8',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-False_stop-False.svg': '89399497a7a9ccf4a3bed6c5e9f62a28',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-False_stop-True.svg': '188d09979d348f1c344cf4f93d919ca1',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-True_stop-False.svg': 'fcd9fdd39d8d6fcd79fc937bf60a48b1',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-False_gtoa-True_stop-True.svg': '195d6bfc06d6db66fc1a9088e8d3ae0a',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-False_stop-False.svg': 'b6171a7e62b3a3e8df76f47aa12b23fa',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-False_stop-True.svg': '24e9514794c9c9edada080f79730d5f0',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-True_stop-False.svg': '82370dba89998bc84c2bf135c44062f5',
            'mutation/Tests/Mutation/zz_nt_title-False_apobec-True_gtoa-True_stop-True.svg': '56a88474ce55bc9827186d7ededc8cda',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-False_stop-False.svg': 'ff3f8bd99185349740f0008ab134451f',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-False_stop-True.svg': '5458d42538c1d7a5e4eb099c021f9ccc',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-True_stop-False.svg': 'd3aa4c4244668c5d126d69f0f1574cd9',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-False_gtoa-True_stop-True.svg': '54fa2ee45bbde283efd5682cf2c652e3',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-False_stop-False.svg': 'b5956a5593fbcefa702fe135b98cfafe',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-False_stop-True.svg': 'e31a381eebb6d839f86411f89cf40dad',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-True_stop-False.svg': 'fe8397e108a06237339ba9a776115fe0',
            'mutation/Tests/Mutation/zz_nt_title-True_apobec-True_gtoa-True_stop-True.svg': '95bd4713eae9aa6cb90b51c2399ea65f'
        }

        #new_hashes: dict = {}

        # PDFs have a creation date in them, so can't be checked by hash as they're different each time
        formats = ["PS", "EPS", "SVG", "JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"] # , "PDF"

        for format in formats:
            if not self.extended_formats and format in ["JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"]:
                print(f"Skipping {format}")
                continue

            with self.subTest(format=format):
                file_name: str = f"mutation/Tests/Mutation/zz_nt.{format.lower()}"
                self.mutation_plot_nt.draw_mismatches(file_name, output_format=format)
                #new_hashes[file_name] = file_hash(file_name=file_name)
                #print(f"\n{file_name}")
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
                #new_hashes[file_name] = file_hash(file_name=file_name)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        #print(f"Extended: {self.extended_formats}")
        #pprint(new_hashes)

    def test_mutation_plot_creates_valid_plot_aa(self):
        """ MutationPlot should create a valid plot """

        # return

        hashes: dict = {
            'mutation/Tests/Mutation/zz_aa.bmp': 'd038e467802cbee4cf57c7c760859626',
            'mutation/Tests/Mutation/zz_aa.eps': '6f8c9cfd8d68773c8d0e7fda3bf0cf0f',
            'mutation/Tests/Mutation/zz_aa.gif': '245d88d8ce2c1a3bb23b82da1e24cbe2',
            'mutation/Tests/Mutation/zz_aa.jpg': 'cb545a04b4e0221bc8ad1179403432cf',
            'mutation/Tests/Mutation/zz_aa.png': 'b5d55844b1a998198fe5e7ce04169615',
            'mutation/Tests/Mutation/zz_aa.ps': '6f8c9cfd8d68773c8d0e7fda3bf0cf0f',
            'mutation/Tests/Mutation/zz_aa.svg': 'b34d86ba82ea5bfcc0ff3a62d16ad0fe',
            'mutation/Tests/Mutation/zz_aa.tif': '448e0e66e011cf7ac78cefcdfb668db9',
            'mutation/Tests/Mutation/zz_aa.tiff': 'bee2bcacc14699e4e105fa831c009fe5',
            'mutation/Tests/Mutation/zz_aa_title-False_glycosylation-False.svg': 'b34d86ba82ea5bfcc0ff3a62d16ad0fe',
            'mutation/Tests/Mutation/zz_aa_title-False_glycosylation-True.svg': '94e1e63ae7b90b95eb10a7ffa2de3976',
            'mutation/Tests/Mutation/zz_aa_title-True_glycosylation-False.svg': '271c0739a15cd925bfa67253ca9a380c',
            'mutation/Tests/Mutation/zz_aa_title-True_glycosylation-True.svg': 'd7a4377fe71868b905c39100d9dca91a'
        }

        #new_hashes: dict = {}

        # PDFs have a creation date in them, so can't be checked by hash as they're different each time
        formats = ["PS", "EPS", "SVG", "JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"] # , "PDF"

        for format in formats:
            if not self.extended_formats and format in ["JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"]:
                continue

            with self.subTest(format=format):
                file_name: str = f"mutation/Tests/Mutation/zz_aa.{format.lower()}"
                self.mutation_plot_aa.draw_mismatches(file_name, output_format=format)
                #new_hashes[file_name] = file_hash(file_name=file_name)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        titles: list = [None, "A Title"]
        glycosylations: list = [True, False]

        for title, glycosylation in product(titles, glycosylations):
            with self.subTest(title=title, glycosylation=glycosylation):
                file_name: str = f"mutation/Tests/Mutation/zz_aa_title-{bool(title)}_glycosylation-{glycosylation}.svg"
                self.mutation_plot_aa.draw_mismatches(file_name, title=title, glycosylation=glycosylation)
                #new_hashes[file_name] = file_hash(file_name=file_name)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        #pprint(new_hashes)

if __name__ == '__main__':
    unittest.main()