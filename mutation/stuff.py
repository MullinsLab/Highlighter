import mutations

import Bio
from Bio import AlignIO, Phylo
from Bio.Align import AlignInfo
from Bio.Phylo import NexusIO, BaseTree
from Bio.Graphics import MutationPlot

# align = AlignIO.read('mutation/Tests/Mutation/highlighter.fasta', 'fasta')
# # align = AlignIO.read('mutation/Tests/Mutation/short_test.fasta', 'fasta')
# mutations = AlignInfo.Mutations(align)

# mutation_plot = MutationPlot(align, title="Mismatches compared to reference")
# mutation_plot.draw("test.svg", apobec=True, g_to_a=True, narrow_markers=True)

# exit()

# V704
align = AlignIO.read("mutation/Tests/Private/thing.fasta", "fasta")
tree = Phylo.read("mutation/Tests/Private/thing.nexus.tre", "nexus")

mutation_plot = MutationPlot(align, tree=tree, ruler_major_ticks=7)
mutation_plot.draw("test.svg", apobec=True, g_to_a=True, narrow_markers=True, sort="tree")

# print(mutations.list_mutations(apobec=False, g_to_a=False))
# print(mutations.list_mutations(apobec=True, g_to_a=False))
# print(mutations.list_mutations(apobec=False, g_to_a=True))
# print(mutations.list_mutations(apobec=True, g_to_a=True))

# print(mutations.list_mutations(reference=1, apobec=False, g_to_a=False))
# print(mutations.list_mutations(reference=1, apobec=True, g_to_a=False))
# print(mutations.list_mutations(reference=1, apobec=False, g_to_a=True))
# print(mutations.list_mutations(reference=1, apobec=True, g_to_a=True))



print("done")