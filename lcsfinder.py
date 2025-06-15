from sequencealignment import *
import pprint

class LCSFinder:
    def __init__(self, seq1, seq2, seq3 = None):
        self.seq1 = seq1
        self.seq2 = seq2
        self.seq3 = seq3
        
        self.lcs_seq = "" # Inicializa o LCS para esta instância
        self.flag = False # Inicializa a flag para esta instância
        # Flag serve para ver se as sequencias entraram num dos computes

        if seq3 == None:
            self.sequences = [seq1, seq2]
        else:   
            self.sequences = [seq1, seq2, seq3]

    #----------> needleman wunsch
    def needleman_wunsch_2(self, seqa, seqb):
        n = seqa.length()
        m = seqb.length()
        match_score = 10
        missmatch_score = -2
        gap_penalty = -5

        # Criar matriz de zeros de dimensoes (n+1) x (m+1) para os scores e para as direcoes
        matrix_score = [[0 for j in range(m+1)] for i in range(n+1)]
        matrix_direction = [[[] for j in range(m+1)] for i in range(n+1)]

        # i anda pelas linhas, j anda pelas colunas
        for i in range(0,n+1):
            for j in range(0,m+1):
                matrix_score[0][j] = gap_penalty * j  # Preencher primeira linha com gaps
                matrix_score[i][0] = gap_penalty* i  # Preencher primeira coluna com gaps:

                matrix_direction[0][j] = ["l"] # Preencher primeira linha com lefts
                matrix_direction[i][0] = ["u"] # Preencher primeira coluna com ups
        
        # Preencher matriz scores : Match = 10, Mismatch = -2 e Gap = -5
        for i in range(1,n+1):
            for j in range(1,m+1):    
                if seqa.char_at(i-1) == seqb.char_at(j-1):
                    score = match_score
                else:
                    score = missmatch_score   
                matrix_score[i][j] = max((matrix_score[i-1][j] + gap_penalty), (matrix_score[i][j-1] + gap_penalty), (matrix_score[i-1][j-1] + score))
                
                # Preencher matriz direcoes
                if (matrix_score[i-1][j] + gap_penalty) >= (matrix_score[i][j-1] + gap_penalty) and (matrix_score[i-1][j] + gap_penalty) >= (matrix_score[i-1][j-1] + score):
                    matrix_direction[i][j] = ["u"]
                if (matrix_score[i][j-1] + gap_penalty) >= (matrix_score[i-1][j] + gap_penalty) and (matrix_score[i][j-1] + gap_penalty) >= (matrix_score[i-1][j-1] + score):
                    matrix_direction[i][j] = ["l"]
                if (matrix_score[i-1][j-1] + score) >= (matrix_score[i][j-1] + gap_penalty) and (matrix_score[i-1][j-1] + score) >= (matrix_score[i-1][j] + gap_penalty):
                    matrix_direction[i][j] = ["d"]

        
        aligned_seqa, aligned_seqb = self.recursive_seq_align_2(seqa, seqb, matrix_direction,(n+1)-1,(m+1)-1,"", "")
        return aligned_seqa, aligned_seqb 
    
    # Função recursiva que descobre o alinhamento com base na matriz        
    
    def recursive_seq_align_2(self,seqa, seqb, matrix_direction,i,j,aligned_seqa,aligned_seqb):
        # Se i ou j forem <= 0, entao chegamos ao fim da matriz, e retornamos as sequencias alinhadas invertidas (comecamos no fim)
        if i <= 0 or j <= 0:
            return aligned_seqa[::-1], aligned_seqb[::-1]
        
        # prioridade: diagonal, depois up, depois left
        if "d" in matrix_direction[i][j]:
            aligned_seqa = aligned_seqa + seqa.char_at(i-1)
            aligned_seqb = aligned_seqb + seqb.char_at(j-1)
            return self.recursive_seq_align_2(seqa, seqb, matrix_direction,i-1,j-1,aligned_seqa,aligned_seqb)
                    
        elif "u" in matrix_direction[i][j]:
            aligned_seqa = aligned_seqa + seqa.char_at(i-1)
            aligned_seqb = aligned_seqb + "-"
            return self.recursive_seq_align_2(seqa, seqb, matrix_direction,i-1,j,aligned_seqa,aligned_seqb)
                
        else: # "l" in matrix_direction[i][j]
            aligned_seqa = aligned_seqa + "-"
            aligned_seqb = aligned_seqb + seqb.char_at(j-1)
            return self.recursive_seq_align_2(seqa, seqb, matrix_direction,i,j-1,aligned_seqa,aligned_seqb)
            
    #<---------- needleman wunsch fim



    def needleman_wunsch_3(self, seqa, seqb, seqc):
        n, m, p = seqa.length(), seqb.length(), seqc.length()
        match_score = 10
        mismatch_score = -2
        gap_penalty = -5

        # Inicializar matriz de pontuação e direção em 3D
        score_matrix = [[[0 for _ in range(p+1)] for _ in range(m+1)] for _ in range(n+1)]
        direction_matrix = [[[[] for _ in range(p+1)] for _ in range(m+1)] for _ in range(n+1)]

        # Inicialização da matriz (fronteiras com gaps)
        for i in range(n+1):
            for j in range(m+1):
                for k in range(p+1):
                    if i == 0 and j == 0 and k == 0:
                        continue
                    scores = []
                    directions = []

                    if i > 0:
                        scores.append(score_matrix[i-1][j][k] + 2 * gap_penalty)
                        directions.append(("u",))
                    if j > 0:
                        scores.append(score_matrix[i][j-1][k] + 2 * gap_penalty)
                        directions.append(("l",))
                    if k > 0:
                        scores.append(score_matrix[i][j][k-1] + 2 * gap_penalty)
                        directions.append(("f",))
                    if i > 0 and j > 0:
                        scores.append(score_matrix[i-1][j-1][k] + gap_penalty + mismatch_score)
                        directions.append(("u", "l"))
                    if i > 0 and k > 0:
                        scores.append(score_matrix[i-1][j][k-1] + gap_penalty + mismatch_score)
                        directions.append(("u", "f"))
                    if j > 0 and k > 0:
                        scores.append(score_matrix[i][j-1][k-1] + gap_penalty + mismatch_score)
                        directions.append(("l", "f"))
                    if i > 0 and j > 0 and k > 0:
                        a, b, c = seqa.char_at(i-1), seqb.char_at(j-1), seqc.char_at(k-1)
                        if a == b == c:
                            score = match_score
                        elif a == b or b == c or a == c:
                            score = mismatch_score  # 2 iguais
                        else:
                            score = 2 * mismatch_score  # todos diferentes
                        scores.append(score_matrix[i-1][j-1][k-1] + score)
                        directions.append(("d",))

                    max_score = max(scores)
                    best_dir = directions[scores.index(max_score)]
                    score_matrix[i][j][k] = max_score
                    direction_matrix[i][j][k] = best_dir

        aligned_seqa, aligned_seqb, aligned_seqc = self.recursive_seq_align_3(seqa, seqb, seqc, direction_matrix, n, m, p, "", "", "")
        return aligned_seqa, aligned_seqb, aligned_seqc


    def recursive_seq_align_3(self, seqa, seqb, seqc, matrix_direction, i, j, k, aligned_a, aligned_b, aligned_c):
        if i <= 0 and j <= 0 and k <= 0:
            return aligned_a[::-1], aligned_b[::-1], aligned_c[::-1]

        direction = matrix_direction[i][j][k][0]

        if direction == "d":
            aligned_a += seqa.char_at(i-1)
            aligned_b += seqb.char_at(j-1)
            aligned_c += seqc.char_at(k-1)
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i-1, j-1, k-1, aligned_a, aligned_b, aligned_c)

        elif direction == "ul":
            aligned_a += seqa.char_at(i-1)
            aligned_b += seqb.char_at(j-1)
            aligned_c += "-"
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i-1, j-1, k, aligned_a, aligned_b, aligned_c)

        elif direction == "uf":
            aligned_a += seqa.char_at(i-1)
            aligned_b += "-"
            aligned_c += seqc.char_at(k-1)
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i-1, j, k-1, aligned_a, aligned_b, aligned_c)

        elif direction == "lf":
            aligned_a += "-"
            aligned_b += seqb.char_at(j-1)
            aligned_c += seqc.char_at(k-1)
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i, j-1, k-1, aligned_a, aligned_b, aligned_c)

        elif direction == "u":
            aligned_a += seqa.char_at(i-1)
            aligned_b += "-"
            aligned_c += "-"
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i-1, j, k, aligned_a, aligned_b, aligned_c)

        elif direction == "l":
            aligned_a += "-"
            aligned_b += seqb.char_at(j-1)
            aligned_c += "-"
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i, j-1, k, aligned_a, aligned_b, aligned_c)

        else:  # direction == "f"
            aligned_a += "-"
            aligned_b += "-"
            aligned_c += seqc.char_at(k-1)
            return self.recursive_seq_align_3(seqa, seqb, seqc, matrix_direction, i, j, k-1, aligned_a, aligned_b, aligned_c)
        #<---------- needleman wunsch 2 fim



    def compute_lcs(self):
        seqs_num = len(self.sequences)
        
        ## Se só existir uma sequência, então não há LCS a calcular
        if seqs_num <= 1:
            print("Error: Not enough sequences provided.")
            return None
        
        ## Se existirem duas ou três sequências, então podemos calcular o LCS
        elif seqs_num == 2:
            self.lcs_seq = self.compute_lcs_2()
            print("LCS sequence: ", self.lcs_seq)

            aligned_seq1, aligned_seq2 = self.needleman_wunsch_2(self.seq1, self.seq2)
            self.flag = True
            return SequenceAlignment(self.seq1, self.seq2, aligned_seq1, aligned_seq2, self.get_lcs_length()) 

        elif seqs_num == 3:
            self.lcs_seq = self.compute_lcs_3()
            self.flag = True
            aligned_seq1, aligned_seq2, aligned_seq3= self.needleman_wunsch_3(self.seq1, self.seq2, self.seq3)
            #needleman recebe objetos Sequence, por isso é que criei estes temporarios
            return SequenceAlignment(self.seq1, self.seq2, aligned_seq1, aligned_seq2, self.get_lcs_length(), self.seq3, aligned_seq3) 
        
        else:
            print("Error: Could not compute LCS algorithm")
            return None
        
    # Funcao da length do LCS
    def get_lcs_length(self):
        if self.flag == False:
            print("LCS not computed yet, therefore no length available.")
            return -1
        else:
            return len(self.lcs_seq)
    
    ### Funções internas

    #----------> Caso de num_seq == 2
    def compute_lcs_2(self):
        lcs_seq = ""

        # como aqui n passei o seq 1 e 2 como argumentos, tenho de usar self.seq1 e self.seq2
        n = self.seq1.length()
        m = self.seq2.length()

        self.n = n
        self.m = m
        
        # Criar matriz de zeros de dimensoes (n+1) x (m+1)
        matrix = [[0 for j in range(m+1)] for i in range(n+1)]
        
        # Preencher matriz
        for i in range(1,n+1):
            for j in range(1,m+1):
                if self.seq1.char_at(i-1) == self.seq2.char_at(j-1):
                    matrix[i][j] = matrix[i-1][j-1] + 1
                else:
                    matrix[i][j] = max(matrix[i-1][j],matrix[i][j-1])
        
        # Demos nome ao resultado do recursive finder, e invertemos a sequencia, pq estava ao contrário
        lcs_seq = self.recursive_finder_2(matrix,(n+1)-1,(m+1)-1,lcs_seq)
        return lcs_seq[::-1]

    # Função recursiva que descobre a maior sequencia com base na matriz
    def recursive_finder_2(self,matrix,i,j,lcs_seq):
        if i <= 0 or j <= 0:
            return lcs_seq

        if self.seq1.char_at(i-1) == self.seq2.char_at(j-1):
            lcs_seq = lcs_seq + self.seq1.char_at(i-1)
            return self.recursive_finder_2(matrix,i-1,j-1,lcs_seq)
        else:
            if matrix[i-1][j] >= matrix[i][j-1]:
                return self.recursive_finder_2(matrix,i-1,j,lcs_seq)
            else: 
                return self.recursive_finder_2(matrix,i,j-1,lcs_seq) 
    #<---------- Caso de num_seq == 2

    
    #----------> Caso de num_seq == 3
    def compute_lcs_3(self):
        lcs_seq = ""

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
        
        # # Imprime matrix preenchida
        # print("LCS Filled Matrix:")
        # for l in matrix:
        #     print(l)
        # print("")
        
        # Descobrir maior sequencia com base na matriz
        lcs_seq = self.recursive_finder_3(matrix,(n+1)-1,(m+1)-1,(r+1)-1,lcs_seq)

        return lcs_seq[::-1]

    def recursive_finder_3(self,matrix,i,j,k,lcs_seq):
        if i <= 0 or j <= 0 or k <= 0:
            return lcs_seq

        if self.seq1.char_at(i-1) == self.seq2.char_at(j-1) and self.seq1.char_at(i-1) == self.seq3.char_at(k-1):
            lcs_seq = lcs_seq + self.seq1.char_at(i-1)
            return self.recursive_finder_3(matrix,i-1,j-1,k-1,lcs_seq)
        else:
            if matrix[i-1][j][k] >= matrix[i][j][k-1] and matrix[i-1][j][k] >= matrix[i][j-1][k]:
                return self.recursive_finder_3(matrix,i-1,j,k,lcs_seq)
            elif matrix[i][j-1][k] >= matrix[i][j][k-1] and matrix[i][j-1][k] >= matrix[i-1][j][k]:
                return self.recursive_finder_3(matrix,i,j-1,k,lcs_seq)
            else: 
                return self.recursive_finder_3(matrix,i,j,k-1,lcs_seq) 
    #<---------- Caso de num_seq == 3