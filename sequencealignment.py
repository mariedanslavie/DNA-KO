class SequenceAlignment:
    seq1 = []
    seq2 = []
    seq3 = []
    sequences = [seq1, seq2] if seq3 == [] else [seq1, seq2, seq3]

    aligned_seq1 = []
    aligned_seq2 = []
    aligned_seq3 = [] 
    aligned_seqs = [aligned_seq1, aligned_seq2] if seq3 == [] else [saligned_seq1, aligned_seq2, aligned_seq3]
    aligned_seq = [0]
    score = 0

    def __init__(self,sequences, aligned_seqs,score):
        self.seq1 = sequences[0]
        self.seq2 = sequences[1]
        self.seq3 = sequences[2] if len(sequences)==3 else None
        self.sequences = sequences 
        self.aligned_seq1 = aligned_seqs[0]
        ##ISTO ESTA MALLL
        self.aligned_seq2 = aligned_seqs[0]
        self.aligned_seq3 = aligned_seqs[2] if len(sequences)==3 else None
        self.aligned_seqs = aligned_seqs
        self.score = score

    # TODO
    def identity(self):
        return self.score
    
    # TODO
    def __str__(self):
        return self.aligned_seq[0]
    
    def create_aligned_seq(self, seq1, seq2, lcs_seq):
        i = j = k = 0
        aligned_seq1 = ""
        aligned_seq2 = ""

        while k < len(lcs_seq):
            # Avança na X até encontrar o próximo caracter da LCS
            while i < len(self.seq1) and self.seq1[i] != lcs_seq[k]:
                aligned_seq1 += self.seq1[i]
                aligned_seq2 += "-"
                i += 1

            # Avança na Y até encontrar o mesmo caracter da LCS
            while j < len(self.seq2) and self.seq2[j] != lcs_seq[k]:
                aligned_seq1 += "-"
                aligned_seq2 += self.seq2[j]
                j += 1

            # Ambos têm o caracter da LCS → alinhamento válido
            if i < len(self.seq1) and j < len(self.seq2) and self.seq1[i] == self.seq2[j] == lcs_seqs[k]:
                aligned_seq1 += self.seq1[i]
                aligned_seq2 += self.seq2[j]
                i += 1
                j += 1
                k += 1

        # Adiciona o que resta de X
        while i < len(self.seq1):
            aligned_seq1 += self.seq1[i]
            aligned_seq2 += "-"
            i += 1

        # Adiciona o que resta de Y
        while j < len(self.seq2):
            aligned_seq1 += "-"
            aligned_seq2 += self.seq2[j]
            j += 1

        return aligned_seq1, aligned_seq2
