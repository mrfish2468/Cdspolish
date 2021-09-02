import pandas as pd
import sys
from Bio import SeqIO
import time
import more_itertools as mit
import numpy as np

def get_nuc(label):
    if label == 0:
        return 'A'
    elif label == 1:
        return 'T'
    elif label == 2:
        return 'C'
    elif label == 3:
        return 'G'


def polish(result, draft_path, path):
    d_position = []
    df = pd.read_feather(result)
    record = SeqIO.read(draft_path,"fasta")

    deletion = df[df['predict'] == 4].position.values #deletion

    insertion = df[(df['predict'] != 5) & (df['predict'] != 4)] #insertion


    d_temp = [list(group) for group in mit.consecutive_groups(deletion)]  
    for key in d_temp:
        if len(key) <= 6:
            d_position = np.concatenate((d_position, key), axis=0)

    seq = []
    seq.append('>{}_CDSpolished\n'.format(record.id))

    for i in range(len(record)): #insertion/deletion
        if i not in d_position: #deletion -> label = 4 
            if i in insertion['position'].values: #insertion
                match = insertion[insertion['position'] == i]     
                seq.append(record[i])  
                for index, row in match.iterrows():                       
                    seq.append(get_nuc(row.predict))
            else:
                seq.append(record[i])

    name = path + '/polished.fasta'
    polished = open(name,'w')
    polished.write(''.join(seq)) 
    polished.write('\n') 
    return name