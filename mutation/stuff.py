import mutations

import Bio
from Bio import AlignIO, Phylo
from Bio.Align import AlignInfo
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Phylo import NexusIO, BaseTree
from Bio.Graphics import MutationPlot


# align = AlignIO.read("mutation/Tests/Private/V704_0011_240-241_REN_highlighter.fasta", "fasta")
# align_info = SummaryInfo(align)
# print(align_info.gap_consensus(threshold=.5))
# if "X" in align_info.gap_consensus():
#     print("X is in consensus")
#     print(f"There are {align_info.gap_consensus(ambiguous='X').count('X')} X's in the consensus")
# exit()


# highlighter_nt match
if True:

    align = AlignIO.read("mutation/Tests/Mutation/highlighter_nt.fasta", "fasta")

    plot = MutationPlot(align)
    plot.draw_matches("highlighter_nt.svg", references = [0, 1])

# V704_0011_240-241_REN
if False:
    align = AlignIO.read("mutation/Tests/Private/V704_0011_240-241_REN_highlighter.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Private/V704_0011_240-241_REN_phyml_tree.tre_nexus.tre", "nexus")

    mutation_plot = MutationPlot(align, tree=tree, top_margin=12, seq_gap=-0.185*2, seq_name_font_size=16, ruler_font_size=12, plot_width=6*72, bottom_margin=45, right_margin=10) # (46*2)-36
    mutation_plot.draw_mismatches("V704_0011_240-241_REN_highlighter-7.svg", apobec=True, g_to_a=True, glycosylation=True, sort="tree", reference="V704_0011_240_REN_NT_29_1", mark_width=7, scheme="ML")
    mutation_plot.draw_mismatches("V704_0011_240-241_REN_highlighter-1.svg", apobec=True, g_to_a=True, glycosylation=True, sort="tree", reference="V704_0011_240_REN_NT_29_1", mark_width=1, scheme="ML")

# V703_0013_090-091_GP_NT_collapsed_by_timepoint
if False:
    align = AlignIO.read("mutation/Tests/Private/V703_0013_090-091_GP_NT_collapsed_by_timepoint.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Private/V703_0013_090-091_GP_NT_collapsed_by_timepoint.phy_phyml_tree.txt_nexus.tre", "nexus")

    # seq_gap=-0.0185*2

    # mutation_plot = MutationPlot(align, tree=tree, top_margin=12, seq_gap=0.9815*2, seq_name_font_size=14, plot_width=6*72)
    mutation_plot = MutationPlot(align, tree=tree, top_margin=12, seq_gap=-0.185*2, seq_name_font_size=16, ruler_font_size=12, plot_width=6*72, bottom_margin=(46*2)-36)
    mutation_plot.draw_mismatches("V703_0013_090-091_GP_NT_collapsed_by_timepoint.svg", apobec=True, g_to_a=True, sort="tree")

# V702_4390_140_env_NT_collapsed_by_TP.fasta
if False:
    align = AlignIO.read("mutation/Tests/Private/V702_4390_140_env_NT_collapsed_by_TP.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Private/V702_4390_140_env_NT_collapsed_by_TP.phy_phyml_tree.txt_nexus.tre", "nexus")

    # mutation_plot = MutationPlot(align, tree=tree, top_margin=12, seq_gap=0.9815*2, seq_name_font_size=14, plot_width=6*72)
    mutation_plot = MutationPlot(align, tree=tree, top_margin=12, seq_gap=-0.0185*2, seq_name_font_size=16, plot_width=6*72)
    mutation_plot.draw_mismatches("V702_4390_140_env_NT_collapsed_by_TP.svg", apobec=True, g_to_a=True, sort="tree")

    # mutation_plot = MutationPlot(align, tree=tree, title="V702_4390_140_env_NT_collapsed_by_TP Offset 0", codon_offset=0)
    # mutation_plot.draw_mismatches("V702_4390_140_env_NT_collapsed_by_TP.offset_0.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

    # mutation_plot = MutationPlot(align, tree=tree, title="V702_4390_140_env_NT_collapsed_by_TP Offset 1", codon_offset=1)
    # mutation_plot.draw_mismatches("V702_4390_140_env_NT_collapsed_by_TP.offset_1.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

    # mutation_plot = MutationPlot(align, tree=tree, title="V702_4390_140_env_NT_collapsed_by_TP Offset 2", codon_offset=2)
    # mutation_plot.draw_mismatches("V702_4390_140_env_NT_collapsed_by_TP.offset_2.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

# Draw plot for hiv_nt.fasta
if False:
    align = AlignIO.read("mutation/Tests/Mutation/hiv_nt.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Mutation/hiv_nt_nexus.tre", "nexus")

    mutation_plot = MutationPlot(align, tree=tree, top_margin=12, seq_gap=-0.185*2, seq_name_font_size=16, ruler_font_size=12, plot_width=6*72, bottom_margin=45, right_margin=10) # (46*2)-36
    mutation_plot.draw_mismatches("DEMO_highlighter.svg", apobec=True, g_to_a=True, sort="tree")

    # mutation_plot = MutationPlot(align, tree=tree, title="HIV_DEMO Offset 0", codon_offset=0)
    # mutation_plot.draw_mismatches("hiv_nt.offset_0.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

    # mutation_plot = MutationPlot(align, tree=tree, title="HIV_DEMO Offset 1", codon_offset=1)
    # mutation_plot.draw_mismatches("hiv_nt.offset_1.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

    # mutation_plot = MutationPlot(align, tree=tree, title="HIV_DEMO Offset 2", codon_offset=2)
    # mutation_plot.draw_mismatches("hiv_nt.offset_2.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

    exit()

if False:
    print(AlignInfo.Mutations.get_mutations(reference='GNS-S', sequence='GNS-S', type='AA', glycosylation=True))
    exit()

if False:
    align = AlignIO.read('mutation/Tests/Mutation/test.fasta', 'fasta')
    mutations = AlignInfo.Mutations(align, type="NT")
    print(mutations.list_mutations(apobec=True, stop_codons=True))
    print(f"Is stop in result: {'Stop' in str(mutations.list_mutations(stop_codons=True))}")
    exit()

if False:
    align = AlignIO.read('mutation/Tests/Mutation/highlighter_aa.fasta', 'fasta')

    mutation_plot = MutationPlot(align, title="Mismatches compared to reference", type="AA")
    mutation_plot.draw_mismatches("test_aa.svg", apobec=True, glycosylation=False)

    exit()

if False:
    print(AlignInfo.Mutations.get_mutations(reference='MRVMEIRRNYQHL--', sequence='MRAMK-RRNYQHL--', type='AA'))
    exit()

if False:
    align = AlignIO.read('mutation/Tests/Mutation/highlighter.fasta', 'fasta')
    # align = AlignIO.read('mutation/Tests/Mutation/short_test.fasta', 'fasta')
    mutations = AlignInfo.Mutations(align)

    mutation_plot = MutationPlot(align, title="Mismatches compared to reference")
    mutation_plot.draw_mismatches("test.svg", apobec=True, g_to_a=True, narrow_markers=True)

    exit()

if False:
    align = AlignIO.read('mutation/Tests/Mutation/highlighter_nt.fasta', 'fasta')
    mutation_plot = MutationPlot(align)
    mutation_plot.draw_mismatches("highlighter_nt.svg", apobec=True, g_to_a=True, stop_codons=True, narrow_markers=True)

    exit()

# V704
if False:
    align = AlignIO.read("mutation/Tests/Private/thing.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Private/thing.nexus.tre", "nexus")

    mutation_plot = MutationPlot(align, tree=tree, ruler_major_ticks=7, title="V704_0011_240_REN_NT")
    mutation_plot.draw_mismatches("V704.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")




print("done")