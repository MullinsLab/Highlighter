import mutations

import Bio
from Bio import AlignIO, Phylo
from Bio.Align import AlignInfo
from Bio.Phylo import NexusIO, BaseTree
from Bio.Graphics import MutationPlot

# Draw plot for hiv_nt.fasta
if True:
    align = AlignIO.read("mutation/Tests/Mutation/hiv_nt.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Mutation/hiv_nt_nexus.tre", "nexus")

    mutation_plot = MutationPlot(align, tree=tree, title="HIV_DEMO")
    mutation_plot.draw("hiv_nt.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")

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
    mutation_plot.draw("test_aa.svg", apobec=True, glycosylation=False)

    exit()

if False:
    print(AlignInfo.Mutations.get_mutations(reference='MRVMEIRRNYQHL--', sequence='MRAMK-RRNYQHL--', type='AA'))
    exit()

if False:
    align = AlignIO.read('mutation/Tests/Mutation/highlighter.fasta', 'fasta')
    # align = AlignIO.read('mutation/Tests/Mutation/short_test.fasta', 'fasta')
    mutations = AlignInfo.Mutations(align)

    mutation_plot = MutationPlot(align, title="Mismatches compared to reference")
    mutation_plot.draw("test.svg", apobec=True, g_to_a=True, narrow_markers=True)

    exit()

if False:
    align = AlignIO.read('mutation/Tests/Mutation/highlighter_nt.fasta', 'fasta')
    mutation_plot = MutationPlot(align)
    mutation_plot.draw("highlighter_nt.svg", apobec=True, g_to_a=True, stop_codons=True, narrow_markers=True)

    exit()

# V704
if False:
    align = AlignIO.read("mutation/Tests/Private/thing.fasta", "fasta")
    tree = Phylo.read("mutation/Tests/Private/thing.nexus.tre", "nexus")

    mutation_plot = MutationPlot(align, tree=tree, ruler_major_ticks=7)
    mutation_plot.draw("test.svg", apobec=True, g_to_a=True, stop_codons=True, sort="tree")




print("done")