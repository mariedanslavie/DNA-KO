from sequence import *


class SequenceDataBase:
    sequences = []

    def __init__(self):
        self.sequences = []


    def add_sequence(self, sequence):
        """Adds a sequence object to the database.
        :param sequence: Sequence - sequence object to be added
        """
        self.sequences.append(sequence)

    def get_sequence_by_id(self, id):
        """Returns a sequence object by its id.
        :param id: str - id of the sequence to be retrieved
        :return: Sequence - sequence object with the given id or None if not found
        """
        # funcao filter recebe uma funcao lambda que verifica se o id do objeto sequence é igual ao id dado
        seq = next(filter(lambda s: s.id == id, self.sequences), None)
        if seq is None:
            print("Error: No sequence found with id: {0}".format(id))
        return seq


    def load_from_fasta(self, filename):
        """For each line in file, create a sequence object with the sequence in the line
        :param filename: str - path to input file with all sequences
        """

        self.sequences = []
        i = 0
        description = ""
        seq = ""
        flag = False

        with open(filename, 'r') as f:
            for line in f:
                if ">" in line:
                    # se sequences estiver vazio, entao o > em questao esta na primeira linha do fasta
                    if flag == False:
                        description = line.strip()
                        flag = True
                        continue

                    else:
                        self.sequences.append(Sequence(str(i), description, seq))
                        seq = ""
                        i = i +1
                        description = line.strip()
                        continue
                else:
                    seq = seq + line.strip()      
            
            # Ainda temos uma sequência a ser adicionada
            self.sequences.append(Sequence(str(i), description, seq))