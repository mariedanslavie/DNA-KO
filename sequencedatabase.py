class SequenceDataBase:
####ISTO PODE ESTAR TUDO ERRADO EU N PERCEBO NADA DISTO MALTA
    def __init__(self):
        self.sequence = {}

###TODO
    def add_sequence(self, sequence):
        SequenceDataBase.append(sequence)

###TODO
    def get_sequence_by_id(self, id):
        return sequence.id

###TODO  ANTIGO
    def load_from_fasta(self, filename):
        with open(file, 'r') as f:
            for line in f:
                sequences.append(Sequence(str(i),"",line.strip()))
                i = i + 1
        return sequences

###TODO ACTUAL 
    def load_from_fasta(self, filename):
        with open(file, 'r') as f:
            for line in f:
                if line contains == <:
                    continue
                else:
                    sequences.append(Sequence(str(i), "", line.strip()))
            return sequences