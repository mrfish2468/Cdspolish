#encoding=utf-8
import os
from os import listdir
from os.path import isfile, isdir, join
import multiprocessing
import time
import sys
from Bio import SeqIO
from modules import alignment
from modules import alignment2csv
from modules import alignment_train
from modules import alignment2csv_train
from modules import prediction
from modules import polishing
from modules.utils.TextColor import TextColor
from modules.utils.FileManager import FileManager

def main():
    assembly = sys.argv[1]
    model_file = sys.argv[2]
    GENUS = sys.argv[3]
    output_dir = sys.argv[4]
    output_dir = output_dir.rstrip()
#    db_path = '/big6_disk/shiuanrung107/cdspolish/seqkit_cdsdb/'
    db_path = '/bip7_disk/hsuenying_108/cdspolish/seqkit_cds_untrimmed/'
    truth_path = '/bip7_disk/hsuenying_108/cdspolish/truth/'

    output_dir = FileManager.handle_output_directory(output_dir)
    contig_output_dir_debug = output_dir + '/debug'
    contig_output_dir_debug = FileManager.handle_output_directory(contig_output_dir_debug)
    assembly_name = assembly.rsplit('/',1)[-1]
    assembly_name = assembly_name.split('.')[0]

    out = []

#get db_size to determine rareness
    size = int(os.path.getsize('{db}{genus}-genus.CDS.fasta'.format(db=db_path, genus=GENUS)))
    size_GB =size/(1024*1024*1024)

    for contig in SeqIO.parse(assembly, 'fasta'):
        contig_output_dir = contig_output_dir_debug + '/' + contig.id
        contig_output_dir = FileManager.handle_output_directory(contig_output_dir)
        contig_name = contig_output_dir + '/' + contig.id +'.fasta'
        SeqIO.write(contig, contig_name, "fasta")

#alignment(minimap+process+pileup)
        db_np = alignment.align(contig_name, db_path, GENUS, contig_output_dir)
        truth_np = alignment_train.align(assembly, truth_path, contig_output_dir)
#transfer to dataframe
        dataframe = alignment2csv.tocsv(contig_name, db_np, truth_np, size_GB, contig_output_dir)
#predict
        result = prediction.predict(dataframe, model_file, contig_output_dir)
#polish
        finish = polishing.polish(result, contig_name, contig_output_dir)
        out.append(finish)

    os.system('cat {} > {}/{}_CDSpolished.fasta'.format(' '.join(out), output_dir, assembly_name))

if __name__ == "__main__":
    main()

