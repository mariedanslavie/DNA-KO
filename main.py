#testing
import sys
import os
from sequence import *
from lcsfinder import *
from sequencealignment import *


### readFile

def readFile(file):
    """For each line in file, create a sequence object with the sequence in the line

    :param file: path to input file with all sequences
    :return: list of sequence objects
    """
    sequences = []
    i = 0
    with open(file, 'r') as f:
        for line in f:
            sequences.append(Sequence(str(i),"",line.strip()))
            i = i + 1

    return sequences


### Main function
def main():

    # Read input file and create list of sequence objects
    if len(sys.argv) >= 2 and os.path.exists(sys.argv[1]):
        sequences = readFile(sys.argv[1])
    else:
        print("Error: Input file not provided or non existent\nUsage: python .\\projeto_1 data\\<INPUT_FILE_PATH>")
        return -1
    
    for seq in sequences:
        print(seq) 
        print(seq.id) 

# Redirects to main
if __name__ == "__main__":
    main()
     