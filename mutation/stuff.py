import mutations

from Bio import AlignIO
from Bio.Align import AlignInfo

align = AlignIO.read('mutation/Tests/Mutation/test.fasta', 'fasta')
print(align[0].id)
mutations = AlignInfo.Mutations(align)

print("done")