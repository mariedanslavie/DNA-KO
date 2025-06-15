class Sequence:
    id = ""
    description = ""
    seq = ""

    def __init__(self, id, description, seq):
        self.id = id
        self.description = description
        self.seq = seq

    def __str__(self):
        """ String representation of the sequence object.
        :return: str - formatted string with sequence id, description and sequence
        """
###### esta not sure if good
        return self.seq
    
    def length(self):
        """ Returns the length of the sequence.
        :return: int - length of the sequence
        """
        return len(self.seq)
    
    def char_at(self,position):
        """ Returns the character at the given position in the sequence.
        :param position: int - position in the sequence
        :return: str - character at the given position or None if position is out of bounds
        """
        if position <= len(self.seq)-1 and position >= -len(self.seq):
            return self.seq[position]
        else:
            print("Given position ({0}) exceeded sequence length ({1})".format(position, len(self.seq)))
            return None