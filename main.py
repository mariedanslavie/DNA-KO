#testing
import sys
import os
from sequence import *
from lcsfinder import *
from sequencealignment import *
from sequencedatabase import *

### readFile -  PRECISO DISTO?

# def readFile(file):
#     """For each line in file, create a sequence object with the sequence in the line

#     :param file: path to input file with all sequences
#     :return: list of sequence objects
#     """
# ### BRO PEGA NAS SEQUENCIAS E METE CADA UMA DENTRO DA LISTA SEQUENCES. O PROBLEMA É QUE A TERCEIRA POSICAO NAO PODE SER VAZIA
# ### PORQUE VaZIA CONTA NA MESMA COMO UMA LISTA.
#     sequences = []
#     i = 0
#     with open(file, 'r') as f:
#         for line in f:
#             sequences.append(Sequence(str(i),"",line.strip()))
#             i = i + 1
#     self.add_sequence
#     return sequences


### Main function
def main():

    database = SequenceDataBase()

    # Read input file and create list of sequence objects
    if len(sys.argv) >= 2 and os.path.exists(sys.argv[1]):
        database.load_from_fasta(sys.argv[1])
    else:
        print("Error: Input file not provided or non existent\nUsage: python .\\projeto_1 data\\<INPUT_FILE_PATH>")
        return -1


    seq1 = database.sequences[0]
    seq2 = database.sequences[1]
    if len(database.sequences) > 2:
        seq3 = database.sequences[2]
    else:
        seq3 = None

    print(seq1)
    print(seq2)
    print(seq3)

 ####### código alternativo #####
 # print do lcs   
    
    if seq3 == None:
        finder = LCSFinder(seq1, seq2)      
    else:
        finder = LCSFinder(seq1, seq2, seq3)


    sequence_alignment = finder.compute_lcs()
    print("The largest found sequence was: {0}".format(finder.lcs_seq))

# ######## código original ###########
#     finder = LCSFinder(sequences)
#     sequence_alignment = finder.compute_lcs()
#     print("The largest found sequence was: {0}".format(finder.lcs_seq))

# Redirects to main
if __name__ == "__main__":
    main()
     