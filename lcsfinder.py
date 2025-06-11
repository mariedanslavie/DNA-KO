from sequencealignment import *
import pprint
### tirei o numpy e o itertools (confirma se era mesmo só preciso para o codigo de n sequencias)


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

#### Se as sequencias entrarem num destes if, a flag fica true
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

#### primeiro (0) é o que se vai meter, o j anda no range de 0 até m e é as posições na matriz -- COLUNAS
#### primeiro (0) é o que se vai meter, o i anda no range de 0 até n e é as posições na matriz -- LINHAS
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

### Matriz bonitinha printed out
        for l in matrix:
            print(l)
        print("")
        
        # Descobrir maior sequencia com base na matriz
        aligned_seq = self.recursive_finder_2(matrix,(n+1)-1,(m+1)-1,aligned_seq)
### Demos nome ao resultado do recursive finder, e invertemos a sequencia, pq estava ao contrário

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