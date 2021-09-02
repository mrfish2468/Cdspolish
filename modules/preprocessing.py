import pandas as pd
import numpy as np

def preprocessing(data):
    new_draft = []
    draft = data['draft']
    for i in draft:
        if i == 'A':
            new_draft.append(0)
        elif i == 'T':
            new_draft.append(1)
        elif i == 'C':
            new_draft.append(2)
        elif i == 'G':
            new_draft.append(3)
        elif i == '-':
            new_draft.append(4)

    data['draft'] = new_draft

    return data
    
def double_ins_chaos(data):
    chaos = np.zeros((data.shape[0],), dtype=int)
    ins_index = data[data.draft == '-'].index.values
    side = 10
    dis = np.diff(ins_index)
    idxs = np.where((dis < 10) & (dis > 1))[0]
    

    for i in idxs:        
        head = data[data.index == ins_index[i]]
        head_total = head[['Ins_A', 'Ins_T', 'Ins_C', 'Ins_G', 'Ins_A_hypo', 'Ins_T_hypo', 'Ins_C_hypo', 'Ins_G_hypo']].sum(axis=1).values[0]
        tail = data[data.index == ins_index[i+1]]
        tail_total = tail[['Ins_A', 'Ins_T', 'Ins_C', 'Ins_G', 'Ins_A_hypo', 'Ins_T_hypo', 'Ins_C_hypo', 'Ins_G_hypo']].sum(axis=1).values[0]
        total_weight = head.total_weight.values[0]
        if (head_total + tail_total) == total_weight:   
            if head_total > tail_total:
                chaos[ins_index[i]] = 1
                chaos[ins_index[i+1]] = 2
            elif head_total < tail_total:
                chaos[ins_index[i]] = 2
                chaos[ins_index[i+1]] = 1
            else:
                chaos[ins_index[i]] = 1
                chaos[ins_index[i+1]] = 1
    data['chaos'] = chaos
    return data