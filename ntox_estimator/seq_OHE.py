import numpy as np

class SeqtoOHE(object):
    def __init__(self):
        pass
        
    def gen_seq_np(self,sequence):

        # one-hot-encoding information
        aa_matrix = {'R':0,'H':1,'K':2,'D':3,'E':4,'S':5,'T':6,'N':7,'Q':8,\
                    'C':9,'G':10,'P':11,'A':12,'V':13,'I':14,'L':15,'M':16,\
                     'F':17,'Y':18,'W':19,'X':20}
        self.seq_id = []
        self.seq = []
        self.np_data = np.zeros((len(sequence),300,21))
        for i,l in enumerate(sequence):
            self.seq_id.append(l.id)
            self.seq.append(str(l.seq))
            for j,alp in enumerate(l.seq):
                try:
                    aa_pos = aa_matrix[alp]
                except:
                    aa_pos=20
                self.np_data[i,j,aa_pos]=1
        return self.np_data, self.seq_id, self.seq