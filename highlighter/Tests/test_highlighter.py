import unittest, os, pathlib

from itertools import product
from pprint import pprint

import highlighter

from Bio import AlignIO, Phylo, SeqUtils
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Graphics import HighlighterPlot

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

class HighlighterStaticMismatchTests(unittest.TestCase):
    """ Test the static mismatch methods of the highlighter object """

    def test_get_mismatches_error_if_no_type_given(self):
        """ get_mismatches should error if no type is provided """

        with self.assertRaises(ValueError):
            AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references=['ATGC'])

    def test_get_mismatches_error_if_no_sequence(self):
        """ get_mismatches should error if no sequence is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_mismatches(references='ATGC', seq_type='NT')

    def test_get_mismatches_error_if_no_reference(self):
        """ get_mismatches should error if no reference is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_mismatches(sequence='ATGC', seq_type='NT')

    def test_get_mismatches_error_if_bad_types(self):
        """ get_mismatches should error if bad types are provided """

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_mismatches(sequence=1, references='ATGC', seq_type='NT')

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references=[], seq_type='NT')

    def test_get_mismatches_error_if_different_lengths(self):
        """ get_mismatches should error if the sequence and reference are different lengths """

        with self.assertRaises(ValueError):
            AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references='ATG', seq_type='NT')

    def test_get_mismatches_succeeds_with_mixed_types(self):
        """ get_mismatches should succeed with mixed types """

        try:
            AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references=Seq('ATGC'), seq_type='NT')
        except:
            self.fail('get_mismatches failed with mixed types')

        try:
            AlignInfo.Highlighter.get_mismatches(sequence=SeqRecord(Seq('ATGC')), references='ATGC', seq_type='NT')
        except:
            self.fail('get_mismatches failed with mixed types')

    def test_get_mismatches_returns_dict(self):
        """ get_mismatches should return a dictionary """

        self.assertEqual(type(AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references='ATGC', seq_type='NT')), dict)

    def test_get_mismatches_returns_empty_dict_given_same_sequence_and_reference(self):
        """ get_mismatches should return an empty dictionary if the sequence and reference are the same """

        self.assertEqual(AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references='ATGC', seq_type='NT'), {})

    def test_get_mismatches_returns_non_empty_dict_given_different_sequence_and_reference(self):
        """ get_mismatches should return a non-empty dictionary if the sequence and reference are different """

        self.assertNotEqual(AlignInfo.Highlighter.get_mismatches(sequence='ATGC', references='ATGG', seq_type='NT'), {})

    def test_get_mismatches_returns_correct_glycosylation_sites(self):
        """ get_mismatches should return a dictionary with the correct glycosylation sites """

        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='GNS-SQ', sequence='GNS-SQ', seq_type='AA', glycosylation=True), {1: ['Glycosylation']})

    def test_get_mismatches_returns_correct_stop_codons(self):
        """ get_mismatches should return a dictionary with the correct stop codons """

        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='GTAA-', sequence='GTAA-', seq_type='NT', stop_codons=True), {})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='GGGTAA-', sequence='GGGTAA-', seq_type='NT', stop_codons=True), {3: ['Stop codon']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='GG--GTAA-', sequence='GG--GTAA-', seq_type='NT', stop_codons=True), {5: ['Stop codon']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='GGGUAG-', sequence='GGGUAG-', seq_type='NT', stop_codons=True), {3: ['Stop codon']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='TAG', sequence='TAG', seq_type='NT', stop_codons=True), {0: ['Stop codon']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='TGA', sequence='TGA', seq_type='NT', stop_codons=True), {0: ['Stop codon']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='TGATAG', sequence='TGATAG', seq_type='NT', stop_codons=True), {0: ['Stop codon'], 3: ['Stop codon']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='TGATAG', sequence='TGATAG', seq_type='NT', stop_codons=False), {})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='TGATAG', sequence='TGATAG', seq_type='NT'), {})

    def test_get_mismatches_returns_correct_dict_given_different_sequence_and_reference_nt(self):
        """ get_mismatches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(g_to_a=True, references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A', 'G->A mutation'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(apobec=True, references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A'], 6: ['Gap'], 7: ['T']})
        self.assertEqual(AlignInfo.Highlighter.get_mismatches(g_to_a=True, apobec=True, references='GTGCGGC-', sequence='AATGCA-T', seq_type='NT'), {0: ['A', 'G->A mutation', 'APOBEC'], 1: ['A'], 2: ['T'], 3: ['G'], 4: ['C'], 5: ['A', 'G->A mutation'], 6: ['Gap'], 7: ['T']})

    def test_get_mismatches_returns_correct_dict_given_different_sequence_and_reference_aa(self):
        """ get_mismatches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Highlighter.get_mismatches(references='MRVMEIRRNYQHL--', sequence='MRAMK-RRNYQHL--', seq_type='AA'), {2: ['A'], 4: ['K'], 5: ['Gap']})


class HighlighterStaticMatchTests(unittest.TestCase):
    """ Test the static match methods of the highlighter object """

    def test_get_matches_error_if_no_type_given(self):
        """ get_matches should error if no type is provided """

        with self.assertRaises(ValueError):
            AlignInfo.Highlighter.get_matches(sequence='ATGC', references='ATGC')

    def test_get_matches_error_if_no_sequence(self):
        """ get_matches should error if no sequence is provided """

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_matches(references='ATGC', seq_type='NT')

    def test_get_matches_error_if_no_reference(self):
        """ get_matches should error if no reference is provided """
    
        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_matches(sequence='ATGC', seq_type='NT')

    def test_get_matches_error_if_bad_types(self):
        """ get_matches should error if bad types are provided """

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_matches(sequence=1, references='ATGC', seq_type='NT')

        with self.assertRaises(TypeError):
            AlignInfo.Highlighter.get_matches(sequence='ATGC', references=1, seq_type='NT')

    def test_get_matches_error_if_different_lengths(self):
        """ get_matches should error if the sequence and reference are different lengths """

        with self.assertRaises(ValueError):
            AlignInfo.Highlighter.get_matches(sequence='ATGC', references='ATG', seq_type='NT')

    def test_get_matches_succeeds_with_mixed_types(self):
        """ get_matches should succeed with mixed types """

        try:
            AlignInfo.Highlighter.get_matches(sequence='ATGC', references=Seq('ATGC'), seq_type='NT')
        except:
            self.fail('get_matches failed with mixed types')

        try:
            AlignInfo.Highlighter.get_matches(sequence=SeqRecord(Seq('ATGC')), references='ATGC', seq_type='NT')
        except:
            self.fail('get_matches failed with mixed types')

    def test_get_matches_returns_dict(self):
        """ get_matches should return a dictionary """

        self.assertEqual(type(AlignInfo.Highlighter.get_matches(sequence='ATGC', references='ATGC', seq_type='NT')), dict)

    def test_get_matches_returns_non_empty_dict_given_different_sequence_and_reference(self):
        """ get_mismatches should return a non-empty dictionary if the sequence and reference are different """

        self.assertNotEqual(AlignInfo.Highlighter.get_matches(sequence='ATGC', references='ATGG', seq_type='NT'), {})

    def test_get_matches_returns_correct_matches_given_different_sequence_and_reference(self):
        """ get_matches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Highlighter.get_matches(references='GTGCGGC-', sequence='GTGCGGCT', seq_type='NT'), {7: ['Unique']})

    def test_get_matches_returns_correct_matches_given_multiple_references(self):
        """ get_matches should return a dictionary with the correct keys and values """

        self.assertEqual(AlignInfo.Highlighter.get_matches(references=['GTGCGGC-', 'GTGTGGCT'], sequence='GTGCCGCT', seq_type='NT'), {3: [0], 4: ['Unique'], 7: [1]})


class HighlighterObjectTests(unittest.TestCase):
    """ Tests that use the highlighter object """

    def setUp(self):
        """ Set up an align object to use for testing """

        self.short_align = AlignIO.read(pathlib.PurePath(pathlib.Path(__file__).parent.resolve(), 'Highlighter/short_test_nt.fasta'), 'fasta')
        self.short_highlighters = AlignInfo.Highlighter(self.short_align, seq_type='NT')

        self.medium_align = AlignIO.read(pathlib.PurePath(pathlib.Path(__file__).parent.resolve(), 'Highlighter/medium_test_nt.fasta'), 'fasta')
        self.medium_highlighters = AlignInfo.Highlighter(self.medium_align, seq_type='NT')

    def test_list_mismatches_raises_error_on_bad_reference(self):
        """ list_mismatches should raise an error if the reference is not in the alignment """

        with self.assertRaises(IndexError):
            self.short_highlighters.list_mismatches(references='bad_test')

        with self.assertRaises(IndexError):
            self.short_highlighters.list_mismatches(references=7)

    def test_list_matches_raises_error_on_bad_reference(self):
        """ list_matches should raise an error if the reference is not in the alignment """

        with self.assertRaises(IndexError):
            self.short_highlighters.list_matches(references='bad_test')

        with self.assertRaises(IndexError):
            self.short_highlighters.list_matches(references=7)
        
        with self.assertRaises(ValueError):
            self.short_highlighters.list_matches(references={})

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
                self.assertEqual(self.short_highlighters.list_mismatches(references=reference, apobec=apobec, g_to_a=g_to_a), result[reference][apobec][g_to_a])

    def test_list_matches_returns_correct_list(self):
        """ list_matches should return a the correct list """

        result: list = [{}, {}, {0: [0], 7: [1]}, {0: [1], 7: [0]}, {0: [1], 3: ['Unique'], 7: [0]}]
        self.assertEqual(self.medium_highlighters.list_matches(references=[0, 1]), result)

        result: list = [{0: [0], 7: [0]}, {0: [1], 7: [1]}, {0: [0], 7: [1]}, {0: [1], 7: [0]}, {0: [1], 3: ['Unique'], 7: [0]}]
        self.assertEqual(self.medium_highlighters.list_matches(references=[Seq("GTGCGGC_TTTT"), Seq("ATGCGGCTTTTT")]), result)


class HighlighterPlotStaticTests(unittest.TestCase):
    """ Tests that use the static methods of the HighlighterPlot class """

    def test_highlighter_plot_guesses_alignment_type_correctly(self):
        """ HighlighterPlot guesses the alignment type correctly """

        self.assertEqual(HighlighterPlot.guess_alignment_type(AlignIO.read('highlighter/Tests/Highlighter/highlighter_nt.fasta', 'fasta')), 'NT')
        self.assertEqual(HighlighterPlot.guess_alignment_type(AlignIO.read('highlighter/Tests/Highlighter/highlighter_aa.fasta', 'fasta')), 'AA')

    def test_highlighter_plot_significant_digits(self):
        """ HighlighterPlot significant_digits returns correctly """

        self.assertEqual(HighlighterPlot.significant_digits(1234), 1200)
        self.assertEqual(HighlighterPlot.significant_digits(51234), 51000)


class HighlighterPlotTests(unittest.TestCase):
    """ Holds test of the HighlighterPlot class """

    def setUp(self):
        """ Set up an align object to use for testing """

        self.align_nt = AlignIO.read('highlighter/Tests/Highlighter/highlighter_nt.fasta', 'fasta')
        self.highlighter_plot_nt = HighlighterPlot(self.align_nt)

        self.align_aa = AlignIO.read('highlighter/Tests/Highlighter/highlighter_aa.fasta', 'fasta')
        self.highlighter_plot_aa = HighlighterPlot(self.align_aa)

        self.align_hiv_nt = AlignIO.read('highlighter/Tests/Highlighter/hiv_nt.fasta', 'fasta')
        self.tree_hiv_nt = Phylo.read('highlighter/Tests/Highlighter/hiv_nt_nexus.tre', 'nexus')
        self.highlighter_plot_hiv_nt = HighlighterPlot(self.align_hiv_nt, tree=self.tree_hiv_nt)

        try:
            from reportlab.graphics import renderPM as test
            self.extended_formats: bool = True
            del test
        except ImportError:
            self.extended_formats: bool = False

    def test_highlighter_plot_inits_correctly(self):
        """ HighlighterPlot should initialize correctly """

        self.assertIs(self.highlighter_plot_nt.alignment, self.align_nt)

    def test_highlighter_plot_fails_init_with_bad_tree(self):
        """ HighlighterPlot should fail to initialize if a bad tree is provided """

        with self.assertRaises(TypeError):
            HighlighterPlot(self.align_nt, tree='bad_tree')

    def test_highlighter_plot_draw_matches_fails_init_with_bad_scheme(self):
        """ HighlighterPlot.draw_matches should fail to initialize if a bad scheme is provided """

        with self.assertRaises(TypeError):
            self.highlighter_plot_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", scheme=1)

        with self.assertRaises(ValueError):
            self.highlighter_plot_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", scheme="bad_scheme")

        with self.assertRaises(ValueError):
            self.highlighter_plot_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", scheme={})

        with self.assertRaises(ValueError):
            self.highlighter_plot_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", scheme={"references": {}, "multiple": "#808080", "unique": "#EFE645"})

        with self.assertRaises(ValueError):
            self.highlighter_plot_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", scheme={"references": {"#808080"}, "multiple": "", "unique": "#EFE645"})

        with self.assertRaises(ValueError):
            self.highlighter_plot_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", scheme={"references": {"#808080"}, "multiple": "#808080", "unique": ""})

    def test_highlighter_plot_draw_matches_creates_valid_plot_with_custom_scheme(self):
        """ HighlighterPlot should create a valid plot with a custom scheme """

        return

        self.highlighter_plot_hiv_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", output_format="PNG", sort="tree", references=[0,5], scheme={"references": ["#FF0000", "#537EFF"], "multiple": "#808080", "unique": "#EFE64"})
        self.assertEqual(file_hash(file_name="highlighter/Tests/Highlighter/matches.png"), "1adb6806b68074378ca69415caa1b62f")
        os.remove("highlighter/Tests/Highlighter/matches.png")

        self.highlighter_plot_hiv_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", output_format="PNG", sort="tree", references=[0,5], scheme={"references": ["#FF0000", "#537EFF"], "multiple": "#808080", "unique": "#EFE64"}, sequence_labels=False)
        self.assertEqual(file_hash(file_name="highlighter/Tests/Highlighter/matches.png"), "0a6a5caf101b79164100b896c5ccd09c")
        os.remove("highlighter/Tests/Highlighter/matches.png")

        with open("highlighter/Tests/Highlighter/hiv_seq_1.txt", "r") as file:
            reference_1: Seq = Seq(file.read())

        with open("highlighter/Tests/Highlighter/hiv_seq_2.txt", "r") as file:
            reference_2: Seq = Seq(file.read())

        self.highlighter_plot_hiv_nt.draw_matches("highlighter/Tests/Highlighter/matches.png", output_format="PNG", sort="tree", references=[reference_1,reference_2], scheme={"references": ["#FF0000", "#537EFF"], "multiple": "#808080", "unique": "#EFE64"})
        self.assertEqual(file_hash(file_name="highlighter/Tests/Highlighter/matches.png"), "bb949268fffdb4eb17e975a2eba73cb3")
        os.remove("highlighter/Tests/Highlighter/matches.png")

    def test_highlighter_plot_creates_valid_plot_nt(self):
        """ HighlighterPlot should create a valid plot """

        return

        hashes: dict = {
            'highlighter/Tests/Highlighter/zz_nt.bmp': '794e595bbbf630f026975fe78c673c5b',
            'highlighter/Tests/Highlighter/zz_nt.eps': '51623140f62bac22468dbf4d368376ed',
            'highlighter/Tests/Highlighter/zz_nt.gif': '7354f9a6ab0b02311f49fb28581104e7',
            'highlighter/Tests/Highlighter/zz_nt.jpg': '23be6b021af8c852abe47c794472ac17',
            'highlighter/Tests/Highlighter/zz_nt.png': 'b924cbfb797bc921f6a819372e6bf828',
            'highlighter/Tests/Highlighter/zz_nt.ps': '51623140f62bac22468dbf4d368376ed',
            'highlighter/Tests/Highlighter/zz_nt.svg': '5b7cfbcb5e8891270047dda6f3905071',
            'highlighter/Tests/Highlighter/zz_nt.tif': '3c1e29491b068280042e28d834bd4c67',
            'highlighter/Tests/Highlighter/zz_nt.tiff': '10a0921d61022e09d4ceefdb0751e7db',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-False_gtoa-False_stop-False.svg': '5b7cfbcb5e8891270047dda6f3905071',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-False_gtoa-False_stop-True.svg': 'a2ee70456e7d3f8364e5b4932c80631c',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-False_gtoa-True_stop-False.svg': '376c0782ba075f3ecdd385fe611ffbfe',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-False_gtoa-True_stop-True.svg': '186f94ad80db4598d760267184f4f0af',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-True_gtoa-False_stop-False.svg': '646cddcefe8dcaee2663d35ab078a262',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-True_gtoa-False_stop-True.svg': '2ee3200e74a80692c6dc012bff83843d',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-True_gtoa-True_stop-False.svg': 'b61dbf2b41a907dad3cd26d1c08abf93',
            'highlighter/Tests/Highlighter/zz_nt_title-False_apobec-True_gtoa-True_stop-True.svg': 'b0d711a4808c9767f937f460f6dd3063',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-False_gtoa-False_stop-False.svg': '09c09b0248483221c41d9b2423efc27c',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-False_gtoa-False_stop-True.svg': '55976f9d4bfca8e900b3b48d6ccd627c',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-False_gtoa-True_stop-False.svg': '2aa1974fffe8188285b185849c0f0f62',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-False_gtoa-True_stop-True.svg': '32b015425739b1bd16c71660e234a861',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-True_gtoa-False_stop-False.svg': 'ba690971da851d95a083e4b0ae1196ad',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-True_gtoa-False_stop-True.svg': '1e104c1d60a79428752b2eed9560e93a',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-True_gtoa-True_stop-False.svg': 'c1e2436cd6c48825958148196468bf7d',
            'highlighter/Tests/Highlighter/zz_nt_title-True_apobec-True_gtoa-True_stop-True.svg': '4583ed3ab37ca31d37f93ae410ab922b'
        }

        #new_hashes: dict = {}

        # PDFs have a creation date in them, so can't be checked by hash as they're different each time
        formats = ["PS", "EPS", "SVG", "JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"] # , "PDF"

        for format in formats:
            if not self.extended_formats and format in ["JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"]:
                print(f"Skipping {format}")
                continue

            with self.subTest(format=format):
                file_name: str = f"highlighter/Tests/Highlighter/zz_nt.{format.lower()}"
                self.highlighter_plot_nt.draw_mismatches(file_name, output_format=format)
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
                file_name: str = f"highlighter/Tests/Highlighter/zz_nt_title-{bool(title)}_apobec-{apobec}_gtoa-{g_to_a}_stop-{stop_codons}.svg"
                self.highlighter_plot_nt.draw_mismatches(file_name, title=title, apobec=apobec, g_to_a=g_to_a, stop_codons=stop_codons)
                #new_hashes[file_name] = file_hash(file_name=file_name)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        #print(f"Extended: {self.extended_formats}")
        #pprint(new_hashes)

    def test_highlighter_plot_creates_valid_plot_aa(self):
        """ HighlighterPlot should create a valid plot """

        return

        hashes: dict = {
            'highlighter/Tests/Highlighter/zz_aa.bmp': '53da1774e6e1650704870b9456438655',
            'highlighter/Tests/Highlighter/zz_aa.eps': '630ca9d5be700e02d656577220d1644e',
            'highlighter/Tests/Highlighter/zz_aa.gif': 'b3d8d8e621a42c58b8bb4d3bc4496be3',
            'highlighter/Tests/Highlighter/zz_aa.jpg': '8df944cbb6206093d752d526e80dcdbf',
            'highlighter/Tests/Highlighter/zz_aa.png': '6ec90f77f3406c67b6c046bd2b9cebe0',
            'highlighter/Tests/Highlighter/zz_aa.ps': '630ca9d5be700e02d656577220d1644e',
            'highlighter/Tests/Highlighter/zz_aa.svg': '74721078610e8ea0abc788be9b933673',
            'highlighter/Tests/Highlighter/zz_aa.tif': '1fe9af0b23b2e3a042d1791da652a90e',
            'highlighter/Tests/Highlighter/zz_aa.tiff': 'fc402e450fd14e6c19a1c108a8e7bf5d',
            'highlighter/Tests/Highlighter/zz_aa_title-False_glycosylation-False.svg': '74721078610e8ea0abc788be9b933673',
            'highlighter/Tests/Highlighter/zz_aa_title-False_glycosylation-True.svg': 'e0b4e7cb74b605dfdefba97c547f5eec',
            'highlighter/Tests/Highlighter/zz_aa_title-True_glycosylation-False.svg': 'e8a9b91b480b6c87cb4d7b44ad0abfe5',
            'highlighter/Tests/Highlighter/zz_aa_title-True_glycosylation-True.svg': '37d21877a6d65ce214d89fd220e21401'
        }

        #new_hashes: dict = {}

        # PDFs have a creation date in them, so can't be checked by hash as they're different each time
        formats = ["PS", "EPS", "SVG", "JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"] # , "PDF"

        for format in formats:
            if not self.extended_formats and format in ["JPG", "BMP", "GIF", "PNG", "TIFF", "TIF"]:
                continue

            with self.subTest(format=format):
                file_name: str = f"highlighter/Tests/Highlighter/zz_aa.{format.lower()}"
                self.highlighter_plot_aa.draw_mismatches(file_name, output_format=format)
                #new_hashes[file_name] = file_hash(file_name=file_name)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        titles: list = [None, "A Title"]
        glycosylations: list = [True, False]

        for title, glycosylation in product(titles, glycosylations):
            with self.subTest(title=title, glycosylation=glycosylation):
                file_name: str = f"highlighter/Tests/Highlighter/zz_aa_title-{bool(title)}_glycosylation-{glycosylation}.svg"
                self.highlighter_plot_aa.draw_mismatches(file_name, title=title, glycosylation=glycosylation)
                #new_hashes[file_name] = file_hash(file_name=file_name)
                self.assertEqual(file_hash(file_name=file_name), hashes[file_name])
                os.remove(file_name)

        #pprint(new_hashes)

if __name__ == '__main__':
    unittest.main()