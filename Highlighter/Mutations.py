from functools import cache

from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Mutations:
    """ Get mutation info from an alignment """

    def __init__(self, alignment):
        """ Initialize the Mutations object """

        self.alignment = alignment
        self.mutations: dict[dict[int: list]] = {}

    def test(self):
        print("test")

    @staticmethod
    def get_mutations(*, sequence: str|Seq|SeqRecord, reference: str|Seq|SeqRecord, apobec: bool=False, g_to_a: bool=False) -> dict[int: list]:
        """ Get mutations from a list of sequences and a reference sequence 
        returns a dictionary of mutations where the key is the position of the mutation and the value is a list of types of mutations """

        if type(sequence) is Seq:
                sequence = str(sequence)
        elif type(sequence) is  SeqRecord:
                sequence = str(sequence.seq)
        elif type(sequence) is not str:
                raise TypeError(f"Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}")

        if type(reference) is Seq:
                reference = str(reference)
        elif type(sequence) is  SeqRecord:
                reference = str(reference.seq)
        elif type(reference) is not str:
                raise TypeError(f"Expected reference to be a string, Seq, or SeqRecord, got {type(reference)}")
            
        if len(sequence) != len(reference):
            raise ValueError("Reference and sequence must be the same length")
        
        return Mutations.get_mutations_cached(sequence=sequence, reference=reference, apobec=apobec, g_to_a=g_to_a)
        

    @cache
    @staticmethod
    def get_mutations_cached(*, sequence: str, reference: str, apobec: bool, g_to_a: bool) -> dict[int: list]:
        """ Get mutations from a sequence and a reference sequence
        separated out so it can be cached (Seq and SeqRecord are not hashable) """

        mutations: dict = {}

        if sequence == reference:
            return mutations

        for base_index in range(len(sequence)):
            if reference[base_index] != sequence[base_index]:
                mutations[base_index] = []
                
                if sequence[base_index] != "-":
                    mutations[base_index].append(sequence[base_index])
                else:
                    mutations[base_index].append("Gap")

                if reference[base_index] == "G" and sequence[base_index] == "A":
                    if g_to_a:
                        mutations[base_index].append("G->A mutation")

                    if apobec and base_index < len(sequence)-3 and sequence[base_index+1] in "AG" and sequence[base_index+2] != "C":
                        mutations[base_index].append("APOBEC")
        
        return mutations
    

AlignInfo.Mutations = Mutations
