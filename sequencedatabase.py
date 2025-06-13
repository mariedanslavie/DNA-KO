from sequence import *


class SequenceDataBase:
    sequences = []

    def __init__(self):
        self.sequences = []

# DEVE ESTAR BOM
    def add_sequence(self, sequence):
        self.sequences.append(sequence)

# DEVE ESTAR BOM
    def get_sequence_by_id(self, id):
        for seq in self.sequence:
            if seq.id == id:
                return seq
            else:
                print("No sequence found with id: {0}".format(id))
                return None

# DEVE ESTAR BOM
    def load_from_fasta(self, filename):
        self.sequences = []
        i = 0
        description = ""
        seq = ""
        flag = False

        with open(filename, 'r') as f:
            for line in f:
                if ">" in line:
    #se sequences estiver vazio, entao o > em questao esta na primeira linha do fasta
                    if flag == False:
                        description = line.strip()
                        flag = True
                        continue

                    else:
                        self.sequences.append(Sequence(str(i), description, seq))
                        i = i +1
                        description = line.strip()
                        
                    ##garantir que nao lemos mais que 3 sequencias
                        if i >= 3:
                            break
                            
                        continue
                else:
                    seq = seq + line.strip()      