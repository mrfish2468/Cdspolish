from Bio import SeqIO
import numpy as np
import pandas as pd
import sys
from collections import Counter

np.seterr(divide='ignore',invalid='ignore')

def tocsv(read, db_np, truth_np, size_GB, path):

    record = SeqIO.read(read,"fasta")
    genome_size = len(record)
    db_arr = np.load(db_np)  
    truth_arr = np.load(truth_np)


    db_arr_hypo, db_arr_nonhypo, db_base, db_Ins_hypo, db_Ins_nonhypo = db_arr['arr_0'], db_arr['arr_1'], db_arr['arr_2'], db_arr['arr_3'], db_arr['arr_4']
    truth_arr, truth_base, truth_Ins = truth_arr['arr_0'], truth_arr['arr_1'], truth_arr['arr_2']

    homo_arr = np.zeros(genome_size, dtype=np.int)
    
    position = []
    draft = []
    A, T, C , G = [], [], [], []
    gap = []
    gap_hypo = []
    ins_A, ins_T, ins_C, ins_G = [], [], [], []
    ins_A_hypo, ins_T_hypo, ins_C_hypo, ins_G_hypo = [], [], [], []
    homo = []
    total_weight = []
    Rareness = []
    label = []

    indel_length=[]

    #identify rareness: common==0 rare==1
    if size_GB > 1:
        rareness = 0
    else:
        rareness = 1
    
    for i in range(genome_size):
        count = 0
        count2=1

        if truth_base[i] != 0 and db_base[i] != 0:
            position.append(i)
            draft.append(record[i])
            A.append((db_arr_hypo[i][0]+db_arr_nonhypo[i][0]))
            T.append((db_arr_hypo[i][1]+db_arr_nonhypo[i][1]))
            C.append((db_arr_hypo[i][2]+db_arr_nonhypo[i][2]))
            G.append((db_arr_hypo[i][3]+db_arr_nonhypo[i][3]))
            gap_hypo.append(db_arr_hypo[i][4])
            gap.append(db_arr_nonhypo[i][4])
            ins_A_hypo.append(0)
            ins_T_hypo.append(0)
            ins_C_hypo.append(0)
            ins_G_hypo.append(0)
            ins_A.append(0)
            ins_T.append(0)
            ins_C.append(0)
            ins_G.append(0)
            total_weight.append(db_base[i])
            Rareness.append(rareness)
            

            index = i
            nuc = record[i]
            while i + 1 < genome_size and record[i+1] == nuc:                               
                count += 1
                i += 1                
            i = index
            
            while record[i] and record[i] == nuc:
                count += 1
                i -= 1
            i = index
            
            homo_arr[i] = count
            
            if db_arr_hypo[i][4]==0 or db_arr_nonhypo[i][4]==0:
                indel_length.append(0)
                
            else:
                index2=i
                while (db_arr_hypo[i][4]+db_arr_nonhypo[i][4])/db_base[i]>0.5 and (db_arr_hypo[i+1][4]+db_arr_nonhypo[i+1][4])/db_base[i+1]>0.5:
                    count2 += 1
                    i += 1
                i=index2
            
                while (db_arr_hypo[i][4]+db_arr_nonhypo[i][4])/db_base[i]>0.5 and (db_arr_hypo[i-1][4]+db_arr_nonhypo[i-1][4])/db_base[i-1]>0.5:
                    count2 += 1
                    i -= 1
                i = index2
                indel_length.append(count2)
                
            
            
#add label_no_correction label_gap
            if truth_arr[i][4] > 0:
                label.append(4)
            else:
                label.append(5)


#check insertion
            for j in range(4):

                db_Ins_A = (db_Ins_hypo[i][j][0]+db_Ins_nonhypo[i][j][0])
                db_Ins_T = (db_Ins_hypo[i][j][1]+db_Ins_nonhypo[i][j][1])
                db_Ins_C = (db_Ins_hypo[i][j][2]+db_Ins_nonhypo[i][j][2])
                db_Ins_G = (db_Ins_hypo[i][j][3]+db_Ins_nonhypo[i][j][3])

                if db_Ins_A != 0 or db_Ins_T != 0 or db_Ins_C != 0 or db_Ins_G != 0:
                    position.append(i)
                    draft.append('-')
                    A.append(0)
                    T.append(0)
                    C.append(0)
                    G.append(0)
                    gap.append(0)
                    gap_hypo.append(0)
                    Rareness.append(rareness)
#calculate indel length
                    if j==0:
                        count3=1
                        if ((db_Ins_hypo[i][j+1][0]+db_Ins_nonhypo[i][j+1][0])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][1]+db_Ins_nonhypo[i][j+1][1])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][2]+db_Ins_nonhypo[i][j+1][2])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][3]+db_Ins_nonhypo[i][j+1][3])/db_base[i]>0.5) :
                            count3+=1
                            if ((db_Ins_hypo[i][j+2][0]+db_Ins_nonhypo[i][j+2][0])>0.5 or (db_Ins_hypo[i][j+2][1]+db_Ins_nonhypo[i][j+2][1])/db_base[i]>0.5 or (db_Ins_hypo[i][j+2][2]+db_Ins_nonhypo[i][j+2][2])/db_base[i]>0.5 or (db_Ins_hypo[i][j+2][3]+db_Ins_nonhypo[i][j+2][3])/db_base[i]>0.5) :
                                count3+=1
                                if ((db_Ins_hypo[i][j+3][0]+db_Ins_nonhypo[i][j+3][0])/db_base[i]>0.5 or (db_Ins_hypo[i][j+3][1]+db_Ins_nonhypo[i][j+3][1])/db_base[i]>0.5 or (db_Ins_hypo[i][j+3][2]+db_Ins_nonhypo[i][j+3][2])/db_base[i]>0.5 or (db_Ins_hypo[i][j+3][3]+db_Ins_nonhypo[i][j+3][3])/db_base[i]>0.5) :
                                    count3+=1
                                    indel_length.append(count3)
                                    
                                else:
                                    indel_length.append(count3)
                                    
                            else:
                                indel_length.append(count3)
                                
                                
                        else:
                            indel_length.append(count3)
                            
                    elif j==1:
                        count3=2
                        if ((db_Ins_hypo[i][j+1][0]+db_Ins_nonhypo[i][j+1][0])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][1]+db_Ins_nonhypo[i][j+1][1])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][2]+db_Ins_nonhypo[i][j+1][2])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][3]+db_Ins_nonhypo[i][j+1][3])/db_base[i]>0.5) :
                            count3+=1
                            if ((db_Ins_hypo[i][j+2][0]+db_Ins_nonhypo[i][j+2][0])>0.5 or (db_Ins_hypo[i][j+2][1]+db_Ins_nonhypo[i][j+2][1])/db_base[i]>0.5 or (db_Ins_hypo[i][j+2][2]+db_Ins_nonhypo[i][j+2][2])/db_base[i]>0.5 or (db_Ins_hypo[i][j+2][3]+db_Ins_nonhypo[i][j+2][3])/db_base[i]>0.5) :
                                count3+=1
                                indel_length.append(count3)
                                
                            else:
                                indel_length.append(count3)
                                
                                
                        else:
                            indel_length.append(count3)
                            
                    elif j==2:
                        count3=3
                        if ((db_Ins_hypo[i][j+1][0]+db_Ins_nonhypo[i][j+1][0])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][1]+db_Ins_nonhypo[i][j+1][1])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][2]+db_Ins_nonhypo[i][j+1][2])/db_base[i]>0.5 or (db_Ins_hypo[i][j+1][3]+db_Ins_nonhypo[i][j+1][3])/db_base[i]>0.5) :
                            count3+=1
                            indel_length.append(count3)
                            
                        else:
                            indel_length.append(count3)
                            
                    elif j==3:
                        count3=4
                        indel_length.append(count3)
                        
                        
#add hypo_ins or nonhypo_ins
                    
                    ins_A_hypo.append(db_Ins_hypo[i][j][0])
                    ins_T_hypo.append(db_Ins_hypo[i][j][1])
                    ins_C_hypo.append(db_Ins_hypo[i][j][2])
                    ins_G_hypo.append(db_Ins_hypo[i][j][3])

                    ins_A.append(db_Ins_nonhypo[i][j][0])
                    ins_T.append(db_Ins_nonhypo[i][j][1])
                    ins_C.append(db_Ins_nonhypo[i][j][2])
                    ins_G.append(db_Ins_nonhypo[i][j][3])
                    total_weight.append(db_base[i])

#add label_ins
                    if truth_Ins[i][j][0] == 0 and truth_Ins[i][j][1] == 0 and truth_Ins[i][j][2] == 0 and truth_Ins[i][j][3] == 0 :
                        label.append(5)              
                    else:
                        for k in range(4):
                            if truth_Ins[i][j][k] != 0:
                                label.append(k)
                                break
    for i in range(genome_size):
        if truth_base[i] != 0 and db_base[i] != 0:
            homo.append(homo_arr[i])
            for j in range(4):
                db_Ins_A = (db_Ins_hypo[i][j][0]+db_Ins_nonhypo[i][j][0])
                db_Ins_T = (db_Ins_hypo[i][j][1]+db_Ins_nonhypo[i][j][1])
                db_Ins_C = (db_Ins_hypo[i][j][2]+db_Ins_nonhypo[i][j][2])
                db_Ins_G = (db_Ins_hypo[i][j][3]+db_Ins_nonhypo[i][j][3])

                if db_Ins_A != 0 or db_Ins_T != 0 or db_Ins_C != 0 or db_Ins_G != 0:
                    if max(db_Ins_A,db_Ins_T,db_Ins_C,db_Ins_G) == db_Ins_A and record[i+1]=='A' :
                        homo.append(homo_arr[i+1])
                    elif max(db_Ins_A,db_Ins_T,db_Ins_C,db_Ins_G) == db_Ins_T and record[i+1]=='T' :
                        homo.append(homo_arr[i+1])
                    elif max(db_Ins_A,db_Ins_T,db_Ins_C,db_Ins_G) == db_Ins_C and record[i+1]=='C' :
                        homo.append(homo_arr[i+1])
                    elif max(db_Ins_A,db_Ins_T,db_Ins_C,db_Ins_G) == db_Ins_G and record[i+1]=='G' :
                        homo.append(homo_arr[i+1])
                    else:
                        homo.append(0)     

    dict = {"position": position,
            "draft": draft,
            "A": A,
            "T": T,
            "C": C,
            "G": G,
            "gap": gap,
            "gap_hypo": gap_hypo,
            "Ins_A": ins_A,
            "Ins_T": ins_T,
            "Ins_C": ins_C,
            "Ins_G": ins_G,
            "Ins_A_hypo": ins_A_hypo,
            "Ins_T_hypo": ins_T_hypo,
            "Ins_C_hypo": ins_C_hypo,
            "Ins_G_hypo": ins_G_hypo,
            "total_weight": total_weight,
            "Rareness": Rareness,
            "homopolymer": homo,
            #"deletion_length":deletion_length,
            #"insertion_length":insertion_length,
            "indel_length":indel_length,
            "label": label
        }


    df_path = path + '/{}_train.feather'.format(record.id)
    csv_path = path + '/{}_alignment_truth.csv'.format(record.id)

    df = pd.DataFrame(dict)

    df.to_csv(csv_path)

    df.to_feather(df_path)
    return df_path
