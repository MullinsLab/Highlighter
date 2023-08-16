import mutations

from Bio import AlignIO
from Bio.Align import AlignInfo

from Bio.Graphics import MutationPlot

align = AlignIO.read('mutation/Tests/Mutation/highlighter.fasta', 'fasta')
# align = AlignIO.read('mutation/Tests/Mutation/short_test.fasta', 'fasta')
mutations = AlignInfo.Mutations(align)
mutation_plot = MutationPlot(align)
mutation_plot.draw("test.svg", apobec=True, g_to_a=True)

# print(mutations.list_mutations(apobec=False, g_to_a=False))
# print(mutations.list_mutations(apobec=True, g_to_a=False))
# print(mutations.list_mutations(apobec=False, g_to_a=True))
# print(mutations.list_mutations(apobec=True, g_to_a=True))

# print(mutations.list_mutations(reference=1, apobec=False, g_to_a=False))
# print(mutations.list_mutations(reference=1, apobec=True, g_to_a=False))
# print(mutations.list_mutations(reference=1, apobec=False, g_to_a=True))
# print(mutations.list_mutations(reference=1, apobec=True, g_to_a=True))



print("done")