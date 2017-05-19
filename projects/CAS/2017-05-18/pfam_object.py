class pfam(object):
    """pfam object to store relevant information for a given pfam hit taken 
    from output .pfd files

    Attributes:
        file
        locus_tag 
        start_loc 
        end_loc 
        strand 
        pfam_num 
        pfam_id 
        sequence
    """

    def __init__(self, line):
        """Return a Customer object whose name is *name* and starting
        balance is *balance*."""
        self.file = line[0]
        self.locus_tag = line[1]
        self.start_loc = line[2]
        self.end_loc = line[3]
        self.strand = line[4]
        self.pfam_num = line[5]
        self.pfam_id = line[6]
        self.sequence = line[7]

    # def withdraw(self, amount):
    #     """Return the balance remaining after withdrawing *amount*
    #     dollars."""
    #     if amount > self.balance:
    #         raise RuntimeError('Amount greater than available balance.')
    #     self.balance -= amount
    #     return self.balance
