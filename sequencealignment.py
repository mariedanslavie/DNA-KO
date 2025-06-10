class SequenceAlignment:
    sequences = []
    aligned_seqs = []
    score = 0

    def __init__(self,sequences,aligned_seqs,score):
        self.sequences = sequences
        self.aligned_seqs = aligned_seqs
        self.score = score

    # TODO
    def identity(self):
        return self.score
    
    # TODO
    def __str__(self):
        return self.aligned_seqs[0]
