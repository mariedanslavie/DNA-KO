from sequencealignment import *
import pprint
### tirei o numpy e o itertools (confirma se era mesmo só preciso para o codigo de n sequencias)


class LCSFinder:
    sequences = []
    seq1 = None
    seq2 = None
    seq3 = None

    sequence_alignment = ""
    flag = False

# Flag serve para ver se as sequencias entraram num dos computes
    def __init__(self, seq1, seq2, seq3 = None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.seq3 = seq3

        if seq3 == None:
            self.sequences = [seq1, seq2]
        else:   
            self.sequences = [seq1, seq2, seq3]

            
    def compute_lcs(self):
        for x in self.sequences:
                print(x)

        seqs_num = len(self.sequences)
        if seqs_num <= 1:
            print("Error: Not enough sequences provided.")
            return None
#################### ALTERAR O NOME DA FUNCAO OU VER SE ESTA BOM
        elif seqs_num == 2:
            self.sequence_alignment = self.compute_lcs_2()
            self.flag = True
            aligned_seq1 = create_aligned_seq(self.seq1, self.sequence_alignment)
            aligned_seq2 = funcao(self.seq2, self.sequence_alignment)

            return SequenceAlignment(self.seq1, self.seq2, aligned_seq1, aligned_seq2, 0) 

        elif seqs_num == 3:
            self.sequence_alignment = self.compute_lcs_3()
            self.flag = True
            aligned_seq1 = funcao(self.seq1, self.sequence_alignment)
            aligned_seq2 = funcao(self.seq2, self.sequence_alignment)
            aligned_seq3 = funcao(self.seq3, self.sequence_alignment)

            return SequenceAlignment(self.seq1, self.seq2, aligned_seq1, aligned_seq2, 0, self.seq3, aligned_seq3) 
        
        else:
            print("Error: Could not compute LCS algorithm")
            return None
        

# DEVE ESTAR BOM
    def get_lcs_length(self):
        if flag == False:
            print("LCS not computed yet, therefore no length available.")

            return -1
        else:
            return len(self.sequence_alignment)
    
    # Funções internas

    #----------> Caso de num_seq == 2
    def compute_lcs_2(self):
        sequence_alignment = ""
### como aqui n passei o seq 1 e 2 como argumentos, tenho de usar self.seq1 e self.seq2

        n = self.seq1.length()
        m = self.seq2.length()

#### primeiro (0) é o que se vai meter, o j anda no range de 0 até m e é as posições na matriz -- COLUNAS
#### primeiro (0) é o que se vai meter, o i anda no range de 0 até n e é as posições na matriz -- LINHAS

        matrix = [[0 for j in range(m+1)] for i in range(n+1)]
        
        # Preencher matriz
        for i in range(1,n+1):
            for j in range(1,m+1):
                if self.seq1.char_at(i-1) == self.seq2.char_at(j-1):
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
        sequence_alignment = self.recursive_finder_2(matrix,(n+1)-1,(m+1)-1,sequence_alignment)
### Demos nome ao resultado do recursive finder, e invertemos a sequencia, pq estava ao contrário

        return sequence_alignment[::-1]

    def recursive_finder_2(self,matrix,i,j,sequence_alignment):
        if i <= 0 or j <= 0:
            return sequence_alignment

        if self.seq1.char_at(i-1) == self.seq2.char_at(j-1):
            sequence_alignment = sequence_alignment + self.seq1.char_at(i-1)
            return self.recursive_finder_2(matrix,i-1,j-1,sequence_alignment)
        else:
            if matrix[i-1][j] >= matrix[i][j-1]:
                return self.recursive_finder_2(matrix,i-1,j,sequence_alignment)
            else: 
                return self.recursive_finder_2(matrix,i,j-1,sequence_alignment) 
    #<---------- Caso de num_seq == 2

    
    #----------> Caso de num_seq == 3
    def compute_lcs_3(self):
        sequence_alignment = ""

        n = self.seq1.length()
        m = self.seq2.length()
        r = self.seq3.length()

        matrix = [[[0 for k in range(r+1)] for j in range(m+1)] for i in range(n+1)]
        
        # Preencher matriz
        for i in range(1,n+1):
            for j in range(1,m+1):
                for k in range(1,r+1):
                    if self.seq1.char_at(i-1) == self.seq2.char_at(j-1) and self.seq1.char_at(i-1) == self.seq3.char_at(k-1):
                        matrix[i][j][k] = matrix[i-1][j-1][k-1] + 1
                    else:
                        matrix[i][j][k] = max(matrix[i-1][j][k],matrix[i][j-1][k],matrix[i][j][k-1])
        
        # Imprime matrix preenchida
        print("LCS Filled Matrix:")
        for l in matrix:
            print(l)
        print("")
        
        # Descobrir maior sequencia com base na matriz
        sequence_alignment = self.recursive_finder_3(matrix,(n+1)-1,(m+1)-1,(r+1)-1,sequence_alignment)

        return sequence_alignment[::-1]

    def recursive_finder_3(self,matrix,i,j,k,sequence_alignment):
        if i <= 0 or j <= 0 or k <= 0:
            return sequence_alignment

        if self.seq1.char_at(i-1) == self.seq2.char_at(j-1) and self.seq1.char_at(i-1) == self.seq3.char_at(k-1):
            sequence_alignment = sequence_alignment + self.seq1.char_at(i-1)
            return self.recursive_finder_3(matrix,i-1,j-1,k-1,sequence_alignment)
        else:
            if matrix[i-1][j][k] >= matrix[i][j][k-1] and matrix[i-1][j][k] >= matrix[i][j-1][k]:
                return self.recursive_finder_3(matrix,i-1,j,k,sequence_alignment)
            elif matrix[i][j-1][k] >= matrix[i][j][k-1] and matrix[i][j-1][k] >= matrix[i-1][j][k]:
                return self.recursive_finder_3(matrix,i,j-1,k,sequence_alignment)
            else: 
                return self.recursive_finder_3(matrix,i,j,k-1,sequence_alignment) 
    #<---------- Caso de num_seq == 3