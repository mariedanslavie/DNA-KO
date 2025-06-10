from sequencealignment import *
import pprint
import numpy as np
import itertools

class LCSFinder:
    sequences = []
    aligned_seq = ""
    flag = False

    def __init__(self,sequences):
        self.sequences = sequences

    def compute_lcs(self):
        seqs_num = len(self.sequences)
        if seqs_num <= 1:
            print("Error: Not enough sequences provided.")
            return None
        elif seqs_num == 2:
            self.aligned_seq = self.compute_lcs_2(self.sequences)
        elif seqs_num == 3:
            self.aligned_seq = self.compute_lcs_3(self.sequences)
        elif seqs_num > 3:
            self.aligned_seq = self.compute_lcs_n(self.sequences)

        flag = True
        return SequenceAlignment(self.sequences, [self.aligned_seq], 0) 
    
    # TODO
    def get_lcs_length(self):
        if self.flag == False:
            return -1
        else:
            return len(self.aligned_seq)
    
    # Funções internas

    #----------> Caso de num_seq == 2
    def compute_lcs_2(self,sequences):
        aligned_seq = ""

        n = sequences[0].length()
        m = sequences[1].length()

        matrix = [[0 for j in range(m+1)] for i in range(n+1)]
        
        # Preencher matriz
        for i in range(1,n+1):
            for j in range(1,m+1):
                if sequences[0].char_at(i-1) == sequences[1].char_at(j-1):
                    matrix[i][j] = matrix[i-1][j-1] + 1
                else:
                    matrix[i][j] = max(matrix[i-1][j],matrix[i][j-1])
        
        # Imprime matrix preenchida
        print("LCS Filled Matrix:")
        for l in matrix:
            print(l)
        print("")
        
        # Descobrir maior sequencia com base na matriz
        aligned_seq = self.recursive_finder_2(matrix,(n+1)-1,(m+1)-1,aligned_seq)

        return aligned_seq[::-1]

    def recursive_finder_2(self,matrix,i,j,aligned_seq):
        if i <= 0 or j <= 0:
            return aligned_seq

        if self.sequences[0].char_at(i-1) == self.sequences[1].char_at(j-1):
            aligned_seq = aligned_seq + self.sequences[0].char_at(i-1)
            return self.recursive_finder_2(matrix,i-1,j-1,aligned_seq)
        else:
            if matrix[i-1][j] >= matrix[i][j-1]:
                return self.recursive_finder_2(matrix,i-1,j,aligned_seq)
            else: 
                return self.recursive_finder_2(matrix,i,j-1,aligned_seq) 
    #<---------- Caso de num_seq == 2

    
    #----------> Caso de num_seq == 3
    def compute_lcs_3(self,sequences):
        aligned_seq = ""

        n = sequences[0].length()
        m = sequences[1].length()
        r = sequences[2].length()

        matrix = [[[0 for k in range(r+1)] for j in range(m+1)] for i in range(n+1)]
        
        # Preencher matriz
        for i in range(1,n+1):
            for j in range(1,m+1):
                for k in range(1,r+1):
                    if sequences[0].char_at(i-1) == sequences[1].char_at(j-1) and sequences[0].char_at(i-1) == sequences[2].char_at(k-1):
                        matrix[i][j][k] = matrix[i-1][j-1][k-1] + 1
                    else:
                        matrix[i][j][k] = max(matrix[i-1][j][k],matrix[i][j-1][k],matrix[i][j][k-1])
        
        # Imprime matrix preenchida
        print("LCS Filled Matrix:")
        for l in matrix:
            print(l)
        print("")
        
        # Descobrir maior sequencia com base na matriz
        aligned_seq = self.recursive_finder_3(matrix,(n+1)-1,(m+1)-1,(r+1)-1,aligned_seq)

        return aligned_seq[::-1]

    def recursive_finder_3(self,matrix,i,j,k,aligned_seq):
        if i <= 0 or j <= 0 or k <= 0:
            return aligned_seq

        if self.sequences[0].char_at(i-1) == self.sequences[1].char_at(j-1) and self.sequences[0].char_at(i-1) == self.sequences[2].char_at(k-1):
            aligned_seq = aligned_seq + self.sequences[0].char_at(i-1)
            return self.recursive_finder_3(matrix,i-1,j-1,k-1,aligned_seq)
        else:
            if matrix[i-1][j][k] >= matrix[i][j][k-1] and matrix[i-1][j][k] >= matrix[i][j-1][k]:
                return self.recursive_finder_3(matrix,i-1,j,k,aligned_seq)
            elif matrix[i][j-1][k] >= matrix[i][j][k-1] and matrix[i][j-1][k] >= matrix[i-1][j][k]:
                return self.recursive_finder_3(matrix,i,j-1,k,aligned_seq)
            else: 
                return self.recursive_finder_3(matrix,i,j,k-1,aligned_seq) 
    #<---------- Caso de num_seq == 3

    def compute_lcs_n(self,sequences):
        aligned_seq = ""

        # Obter lista com dimensão de cada sequência: ex: [3,4,3,3]
        seq_lengths = [seq.length() for seq in sequences]
        # Lista de ranges para cada dimensão: [range(0, 3), range(0, 4), range(0, 3), range(0, 3)]
        iterables = [range(1,length) for length in seq_lengths]
        # Matriz com todas as combinações de iterações possiveis para todas as sequências
        product = list(itertools.product(*iterables))

        # Criar tensor de dimensão igual ao número de sequências
        matrix = np.zeros(seq_lengths, dtype=int)

        # Obter numero de sequencias
        num_seqs = len(sequences)

        for i in range(len(product)):
            list_iterador = list(product[i])
            if all([sequences[0].char_at(product[i][0]) == seq.char_at(product[i][t]-1) for seq, t in zip(sequences,product[i])]):
                # Obter iterador da diagonal anterior
                it_diagonal = tuple([x-1 for x in list_iterador])

                matrix[product[i]] = matrix.item(it_diagonal) + 1
            else:
                it_anteriores = [[x-1 if k == n else x for k,x in enumerate(list_iterador)] for n in range(num_seqs)]

                matrix[product[i]] = max([matrix.item(tuple(it)) for it in it_anteriores])

        for i in range(len(product)):
            print(matrix.item(product[i]))