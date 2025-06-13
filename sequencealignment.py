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

#########ISTO É NECESSRIO?
        self.aligned_seq1 = aligned_seq1 
        self.aligned_seq2 = aligned_seq2 
        self.aligned_seq3 = aligned_seq3 

        self.score = score

        if seq3 == None:
            self.sequences = [seq1, seq2]
        else:
            self.sequences = [seq1, seq2, seq3] 

      
 

        ###################### NOSSO CODIGO align 2 ################## 
    def create_aligned_seq2(self, seq2, sequence_alignment):
        j = 0
        k = 0
        positions2 = 0
        while j < len(self.seq2)-1 and j < len(self.lcs_seq)-1:
            if self.seq2.char_at(j) == self.lcs_seq[k]:
                aligned_seq2 = aligned_seq2 + self.seq2.char_at(j)
                positions2 = positions2 + j
                j+=1
                k+=1
            else: 
                aligned_seq2.append('-')
                j+=1
        return aligned_seq2

##### TODO ----- Feito?
    def identity(self):
        positions1 = 0
        positions2 = 0
        if positions1 != positions2:
            print("Error: Sequences are not aligned properly, cannot compute identity.")
            return -1
        else:
            aligned_seq1 = create_aligned_seq1(self.seq1, self.seq2, self.sequence_alignment)
            aligned_seq2 = create_aligned_seq2(self.seq1, self.seq2, self.sequence_alignment)
            #self.score = len(positions1)/len(aligned_seq1) #############aligned seq é uma string ou uma lista?
            #### VERSAO FUNCIONAL da linha anterior
            self.score = (lambda p, a: len(p) / len(a))(positions1, aligned_seq1)

        return self.score
    
    
# DEVE ESTAR BOM
    def __str__(self):
        if self.seq3 == None:
            return self.aligned_seq1 + "\n" + self.aligned_seq2
        else: 
            return self.aligned_seq1 + "\n" + self.aligned_seq2 + "\n" + self.aligned_seq3



   




       

# ###################### CODIGO DO CHATO TODO MARADO ################## 
#     def create_aligned_seq(self, seq1, seq2, sequence_alignment):
#         i = j = k = 0
#         aligned_seq1 = ""
#         aligned_seq2 = ""

#         while k < len(sequence_alignment):
#             # Avança na X até encontrar o próximo caracter da LCS
#             while i < len(self.seq1) and self.seq1[i] != sequence_alignment[k]:
#                 aligned_seq1 += self.seq1[i]
#                 aligned_seq2 += "-"
#                 i += 1

#             # Avança na Y até encontrar o mesmo caracter da LCS
#             while j < len(self.seq2) and self.seq2[j] != sequence_alignment[k]:
#                 aligned_seq1 += "-"
#                 aligned_seq2 += self.seq2[j]
#                 j += 1

#             # Ambos têm o caracter da LCS → alinhamento válido
#             if i < len(self.seq1) and j < len(self.seq2) and self.seq1[i] == self.seq2[j] == sequence_alignment[k]:
#                 aligned_seq1 += self.seq1[i]
#                 aligned_seq2 += self.seq2[j]
#                 i += 1
#                 j += 1
#                 k += 1

#         # Adiciona o que resta de X
#         while i < len(self.seq1):
#             aligned_seq1 += self.seq1[i]
#             aligned_seq2 += "-"
#             i += 1

#         # Adiciona o que resta de Y
#         while j < len(self.seq2):
#             aligned_seq1 += "-"
#             aligned_seq2 += self.seq2[j]
#             j += 1

#         return aligned_seq1, aligned_seq2
