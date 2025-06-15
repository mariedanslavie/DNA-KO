#testing
import sys
import os
from sequence import *
from lcsfinder import *
from sequencealignment import *
from sequencedatabase import *


### Main function
def main():

    database = SequenceDataBase()

    # Le o file de input se existir
    if len(sys.argv) >= 2 and os.path.exists(sys.argv[1]):
        database.load_from_fasta(sys.argv[1])
    else:
        print("Error: Input file not provided or non existent\nUsage: python .\\projeto_1 data\\<INPUT_FILE_PATH>")
        return -1
    
    # define as sequencias a comparar
    seq1 = database.sequences[0]
    seq2 = database.sequences[1]
    # se existir, define a terceira sequencia
    if len(database.sequences) > 2:
        seq3 = database.sequences[2]
    else:
        seq3 = None

    # Se a terceira sequencia nao existir, entao o finder so recebe duas sequencias
    if seq3 == None:
        finder = LCSFinder(seq1, seq2)
    # Se existir, entao o finder recebe as 3 sequencias    
    else:
        finder = LCSFinder(seq1, seq2, seq3)


    sequence_alignment = finder.compute_lcs()
    print("The largest found sequence was:\t {0}".format(finder.lcs_seq))

    # para fzr prints de teste do needleman wunch
    print("Aligned Sequence 1:\t" + sequence_alignment.aligned_seq1)
    print("Aligned Sequence 2:\t" + sequence_alignment.aligned_seq2)
    if seq3 != None:
        print("Aligned Sequence 3:\t" + sequence_alignment.aligned_seq3)

    print("Identity Score for aligned sequences: {0} %".format(sequence_alignment.identity()))

# Redirects to main
if __name__ == "__main__":
    main()
     