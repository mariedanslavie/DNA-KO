class Sequence:
    id = ""
    description = ""
    seq = ""

    def __init__(self, id, description, seq):
        self.id = id
        self.description = description
        self.seq = seq

    def __str__(self):
        return self.seq
    
    def length(self):
        return len(self.seq)
    
    def char_at(self,position):
        if position <= len(self.seq)-1 and position >= -len(self.seq):
            return self.seq[position]
        else:
            print("Given position ({0}) exceeded sequence length ({1})".format(position, len(self.seq)))
            return None