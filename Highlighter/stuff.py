import Mutations

from Bio import AlignIO
from Bio.Align import AlignInfo

align = AlignIO.read('Highlighter/Tests/Mutations/test.fasta', 'fasta')
mutations = AlignInfo.Mutations(align)

print("done")