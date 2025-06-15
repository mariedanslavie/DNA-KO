from sequence import *
from lcsfinder import *

class SequenceAlignment:
    seq1 = []
    seq2 = []
    seq3 = []

    aligned_seq1 = ""
    aligned_seq2 = ""
    aligned_seq3 = "" 
    
    score = 0

    def __init__(self,seq1, seq2, aligned_seq1, aligned_seq2, score, seq3 = None, aligned_seq3 = None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.seq3 = seq3

        self.aligned_seq1 = aligned_seq1 
        self.aligned_seq2 = aligned_seq2 
        self.aligned_seq3 = aligned_seq3 

        self.score = score

        if seq3 == None:
            self.sequences = [seq1, seq2]
        else:
            self.sequences = [seq1, seq2, seq3] 


    def identity(self):
        """Calculates the identity score of the aligned sequences.
        :return: int - identity score as a percentage
        """
        self.score = len (self.aligned_seq1)
        count = 0
        if self.seq3 == None:
            ### a é do aligned_seq1 e b é do aligned_seq2
            for a,b in zip(self.aligned_seq1, self.aligned_seq2):
                if a == b:
                    count = count + 1
        else:
            for a,b,c in zip(self.aligned_seq1, self.aligned_seq2, self.aligned_seq3):
                if a == b and a == c:
                    count = count + 1
        # x é uma função lambda que calcula a percentagem de identidade
        x = lambda p, a: round((p / a) * 100) if a > 0 else 0
        return x (count, self.score)
    
    
# DEVE ESTAR BOM
    def __str__(self):
        """ String representation of the sequence alignment object.
        :return: str - formatted string with aligned sequences
        """
        if self.seq3 == None:
            return self.aligned_seq1 + "\n" + self.aligned_seq2
        else: 
            return self.aligned_seq1 + "\n" + self.aligned_seq2 + "\n" + self.aligned_seq3

