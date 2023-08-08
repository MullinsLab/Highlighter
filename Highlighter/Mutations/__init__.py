from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Mutations(Align.MultipleSeqAlignment):
    """ A class to store mutations in a multiple sequence alignment """


def get_mutations(*, sequence, reference):
    """ Get mutations from a list of sequences and a reference sequence """

    if type(sequence) is Seq:
            sequence = str(sequence)
    elif type(sequence) is  SeqRecord:
            sequence = str(sequence.seq)
    elif type(sequence) is not str:
            raise TypeError(f'Expected sequence to be a string, Seq, or SeqRecord, got {type(sequence)}')

    if type(reference) is Seq:
            reference = str(reference)
    elif type(sequence) is  SeqRecord:
            reference = str(reference.seq)
    elif type(reference) is not str:
            raise TypeError(f'Expected reference to be a string, Seq, or SeqRecord, got {type(reference)}')
        
    if len(sequence) != len(reference):
        raise ValueError('Reference and sequence must be the same length')
    
    # for base_index in range(0, len(sequence)-1):
    #     if base_
        