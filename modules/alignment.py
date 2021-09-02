from Bio import SeqIO
import numpy as np
import sys
import math
import os

LONG_DELETION_LENGTH = 25
DB_COVERAGE=10
np.set_printoptions(threshold=100)

def process_paf(paf_file, output_file):
    with open(output_file, 'w+') as p:
        with open(paf_file, 'r') as f:
            lines = []
            for line in f:
                line = line.split()      
                Query_name = str(line[0])
                Query_name = Query_name.split("|")
                db_id = str(Query_name[1])
                protein_name = str(Query_name[2])

                Query_length = int(line[1])
                Query_start = int(line[2])
                Query_end = int(line[3])
                strand = str(line[4])
                Target_name = str(line[5])
                Target_start = int(line[7])
                Target_end = int(line[8])
                match_base = int(line[9])
                base = int(line[10])
                cigar = str(line[-1])
                query_cover = float(((Query_end-Query_start)/Query_length)*100)
                identity= float((match_base/base)*100)
                hypo = 0
                if "hypothetical" in protein_name or "Hypothetical" in protein_name:
                    hypo = 1
                if "putative" in protein_name or "Putative" in protein_name:
                    hypo = 1
                if "uncharacterised" in protein_name or "Uncharacterised" in protein_name:
                    hypo = 1
                if "uncharacterized" in protein_name or "Uncharacterized" in protein_name:
                    hypo = 1
                if "unknown" in protein_name or "Unknown" in protein_name:
                    hypo = 1
                temp = [db_id, str(Query_length), str(Query_start), str(Query_end), str(query_cover), strand, Target_name, str(Target_start), str(Target_end), str(match_base), str(base), str(identity), str(hypo), cigar]
                lines.append(temp)
            lines.sort(key = lambda x: (float(x[11]), float(x[4]), int(x[1])), reverse = True)
        for item in lines:
            p.write(" ".join(item) + "\n")

def process_cigar(line, db_dict, arr1, Ins, total, total_weight):
    db_id=line[0]
    identity=float(line[11])
    query_cover=float(line[4])
    query_length=int(line[1])
    cigar = line[13]
    cigar = cigar[5:]
    start_pos = int(line[7])

    flag = 0
    longdel_count = 0
    longdel_status = 0
    rank=0
    for i in cigar: # ex. =ATC*c
        db_dict.setdefault(start_pos,{})[db_id]=identity
        if i == 'A' or i == 'a':  # A:0, T:1, C:2, G:3
            base = 0
        elif i == 'T' or i == 't':
            base = 1
        elif i == 'C' or i == 'c':
            base = 2
        elif i == 'G' or i == 'g':
            base = 3
                        
                    
                        
        if i == '=': #match:1
            flag = 1
        elif i == '*': #mismatch:2
            flag = 2
            mismatch = 0
        elif i == '+': #insertion
            flag = 3
            Ins_pos = 0
            over_ins = []
        elif i == '-': #deletion
            flag = 4
            longdel_status = 0
                        
        elif flag == 1:
            longdel_count = 0     
            rank = list(db_dict[start_pos].keys()).index(db_id)
            if rank <=2 :                       
                arr1[start_pos][base] += DB_COVERAGE/math.pow(2,rank) 
                total_weight[start_pos] += DB_COVERAGE/math.pow(2,rank) 
            else:
                arr1[start_pos][base] += DB_COVERAGE/math.pow(2,3) 
                total_weight[start_pos] += DB_COVERAGE/math.pow(2,3)              
            total[start_pos] += 1
            start_pos += 1 
        elif flag == 2: 
            longdel_count = 0
            if mismatch != 1:
                mismatch += 1
            else:
                rank = list(db_dict[start_pos].keys()).index(db_id)
                if rank <=2 :  
                    arr1[start_pos][base] +=  DB_COVERAGE/math.pow(2,rank) 
                    total_weight[start_pos] += DB_COVERAGE/math.pow(2,rank) 
                else:
                    arr1[start_pos][base] +=  DB_COVERAGE/math.pow(2,3)
                    total_weight[start_pos] += DB_COVERAGE/math.pow(2,3) 
                total[start_pos] += 1 
                start_pos += 1
                            
        elif flag == 3:  #insertion
            #+AAAAA
            #-0123
            #01234
            longdel_count = 0
            if Ins_pos < 4:
                rank = list(db_dict[start_pos-1].keys()).index(db_id)
                if rank <=2 : 
                    Ins[start_pos-1][Ins_pos][base] += DB_COVERAGE/math.pow(2,rank)
                else:
                    Ins[start_pos-1][Ins_pos][base] += DB_COVERAGE/math.pow(2,3)
                Ins_pos += 1
                over_ins.append(base)
            elif Ins_pos == 4:
                for x,y in zip(range(4), over_ins):
                    rank = list(db_dict[start_pos-1].keys()).index(db_id)
                    if rank <=2 : 
                        Ins[start_pos-1][x][y] -= DB_COVERAGE/math.pow(2,rank) 
                    else:
                        Ins[start_pos-1][x][y] -= DB_COVERAGE/math.pow(2,3)
                over_ins = []
                                
                                
        elif flag == 4:    #deletion
                           
            if longdel_status == 0:
                longdel_count += 1
            if longdel_count > LONG_DELETION_LENGTH and longdel_status == 0:
                for i in range(1,LONG_DELETION_LENGTH + 1):
                    rank = list(db_dict[start_pos-i].keys()).index(db_id)
                    if rank <=2 : 
                        arr1[start_pos-i][4] -= DB_COVERAGE/math.pow(2,rank)
                        total_weight[start_pos-i] -= DB_COVERAGE/math.pow(2,rank) 
                    else:
                        arr1[start_pos-i][4] -= DB_COVERAGE/math.pow(2,3) 
                        total_weight[start_pos-i] -= DB_COVERAGE/math.pow(2,3)                             
                    total[start_pos-i] -= 1                                                  
                longdel_status = 1
                longdel_count = 0
            elif longdel_status != 1:
                rank = list(db_dict[start_pos].keys()).index(db_id)
                if rank <=2 : 
                    arr1[start_pos][4] += DB_COVERAGE/math.pow(2,rank) 
                    total_weight[start_pos] += DB_COVERAGE/math.pow(2,rank)   
                else:
                    arr1[start_pos][4] += DB_COVERAGE/math.pow(2,3) 
                    total_weight[start_pos] += DB_COVERAGE/math.pow(2,3)   
                                               
                total[start_pos] += 1
            start_pos+=1
    return arr1, total, total_weight, Ins, db_dict

def aligment(filename, Genome_size,genome_id):
    """[summary]
    
    Arguments:
        filename {paf file name} -- db alignment result
        Genome_size {integer} -- TGS assembly only genome
    
    Returns:
        arr1: genome each base [0-3]: ATCG, [4]: deletion
        total: genome each base total alignment result
        Ins: genome each base [0-3]: 1st Ins, 2nd Ins, 3rd Ins; [0-3]: ATCG
    """
    arr_nonhypo = np.zeros((Genome_size,5), dtype=np.float)
    arr_hypo = np.zeros((Genome_size,5), dtype=np.float)
    total = np.zeros(Genome_size, dtype=np.int)
    Ins_nonhypo = np.zeros((Genome_size, 4, 4), dtype=np.float)
    Ins_hypo = np.zeros((Genome_size, 4, 4), dtype=np.float)

    
    total_weight = np.zeros(Genome_size, dtype=np.float)

    line_count = 0  
    over_ins = []  
    db_dict={}
    
    with open(filename, 'r') as f:      
        for line in f:
            line_count += 1 #paf: the number of lines
            line = line.split()      

            hypo = int(line[12])
            start_pos = int(line[7])
            end_pos = int(line[8])
            middle_pos=int((end_pos+start_pos)/2)
            if total[middle_pos]>=DB_COVERAGE:
                continue 

            if hypo == 1:#hypothetical protein
                arr_hypo, total, total_weight, Ins_hypo, db_dict = process_cigar(line, db_dict, arr_hypo, Ins_hypo, total, total_weight)
            else:#nonhypothetical protein
                arr_nonhypo, total, total_weight, Ins_nonhypo, db_dict = process_cigar(line, db_dict, arr_nonhypo, Ins_nonhypo, total, total_weight)

        iden = '{}_identity.txt'.format(genome_id)
        depth = '{}_depth.txt'.format(genome_id)
        fp=open(iden,'w')   
        for k,v in sorted(db_dict.items()):
              fp.write(str(k)+' '+str(v)+'\n')
        fp.close()

        return arr_hypo, arr_nonhypo, total_weight, Ins_hypo, Ins_nonhypo, #, cds

def align(draft_path, db_path, GENUS, path):
    record = SeqIO.read(draft_path,"fasta")
    Genome_size = len(record)

    paf = '{}/{}.paf'.format(path, record.id)
    npz = '{}/{}.npz'.format(path, record.id)
    processed_paf = '{}/{}_processed.paf'.format(path, record.id)
    genome_id=record.id

#run_minimap
    minimap2_cmd = 'minimap2 -cx asm20 --cs=long -t 32 {draft} {db}{genus}-genus.CDS.fasta> {paf}'\
        .format(draft=draft_path, db=db_path, genus=GENUS, paf=paf)
    os.system(minimap2_cmd)
    if os.stat(paf).st_size == 0: #minimap2 can't align return false
        return False

#process_paf(add idendity and query_cover + sort by idendity query_cover query length)
    process_paf(paf, processed_paf)

#pileup
    arr_hypo, arr_nonhypo, total_weight, Ins_hypo, Ins_nonhypo = aligment(processed_paf, Genome_size, genome_id)
    np.savez(npz, arr_hypo, arr_nonhypo, total_weight, Ins_hypo, Ins_nonhypo)
    
    return npz

