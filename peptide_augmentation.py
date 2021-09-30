from Bio import SeqIO
import numpy as np
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from Bio import SearchIO
import pandas as pd
import re
import math
import argparse
import os
import configparser

config = configparser.ConfigParser()
config.read('aug_config.ini')

blast_bin_path = config['DEFAULT']['BLAST_BIN_PATH']

parser = argparse.ArgumentParser(description='Peptide data augmentation')
parser.add_argument('--fasta',metavar='fasta', type=str, help='Fasta input file')
parser.add_argument('--cpu',metavar='fasta', type=str, help='Max cpu usage (Defalut: 1)', default=1)
parser.add_argument('--eval',metavar='fasta', type=str, help='E-value cutoff (Default: 1E-5)', default=1e-5)

# parser.add_argument('--output',metavar='output', type=str, help='Output file')
args = parser.parse_args()
cpu_usage = args.cpu
evalue = args.eval
def seq_prob_list(alp):
    seq_list = []
    aa_group = [['R','H','K'],\
            ['D','E'],\
            ['S','T','N','Q'],\
            ['C','G','P'],\
            ['A','V','I','L','M','F','Y','W']]
    if random.uniform(0,1)<0.6 and alp!='X':
        for l in aa_group:
            if alp in l:
                l.remove(alp)
                return random.sample(l,1)
    else:
        for l in aa_group:
            if not alp in l:
                seq_list+=l

        return random.sample(seq_list,1)
            


def rand_seq_gen(sequence):

    substitute_count = 0
    miss_count=-1
    seq = []
    substitute_count = 0
    miss_count+=1
    for i, alp in enumerate(sequence):
        if alp!='*':
            ## Random AA addition
            if random.uniform(0,1)<0.1:
                seq.extend(seq_prob_list('X'))
                substitute_count+=1

            ## Random AA substitution
            if i==0 and alp=='M':
                seq.extend(alp)
            elif random.uniform(0,1)<0.5:
                seq.extend(seq_prob_list(alp))
                substitute_count+=1
            else:
                seq.extend(alp)
    
    return ''.join(seq)

def fol_create(fol_name):
    if not os.path.isdir(fol_name):
        os.makedirs(fol_name)


## result file folder create
fol_create('./aug_data')



## BLAST db folder create
fol_create('./aug_data/blast_db')


    
print('Make blast db -------------------------------------------')

                
subprocess.call([blast_bin_path+'/makeblastdb','-dbtype','prot',\
                 '-in',str(args.fasta),\
                 '-input_type','fasta',\
                 '-out','./aug_data/blast_db/blast_db'])
print('Make blast db done --------------------------------------')
print('Read original sequence data file-------------------------')
orig_seq_data = SeqIO.parse(args.fasta,'fasta')
data = []
id_info = []
for l in orig_seq_data:
    data.append(l)
    id_info.append([l.id, len(l.seq), 0])
df_id = pd.DataFrame(id_info)
df_id.columns = ['ID','Seq len','Count']
p = re.compile('([\S]+)#[0-9]+#[0-9]+')


### Parameters for augmentation
aug_pep_num = 10
rand_len = 50
max_trial = 30
blast_cpu_usage = cpu_usage
e_val = evalue
###################################

blast_matched_seq = []
num = -1
res = 0

## randomly generated seq folder create
fol_create('./aug_data/rand_gen_seq')
fol_create('./aug_data/blast_out')




while res<np.shape(df_id)[0] and num<30:
    num += 1
    
    print('Try ',str(num),'/',str(max_trial))
    print('Random sequence generation start ------------------------')
    rand_gen_seq_list = []
    for l in data:
        if (df_id.loc[df_id['ID']==l.id,'Count'].values<=aug_pep_num)[0]:
            for j in range(rand_len):
#                 print(l.seq)
                seq_rec = SeqRecord(Seq(rand_seq_gen(l.seq)), id=l.id+'#'+str(num)+'#'+str(j),\
                                    description='Rand gen seq', name=l.id+'_'+str(j))
                rand_gen_seq_list.append(seq_rec)
    SeqIO.write(rand_gen_seq_list,'./aug_data/rand_gen_seq/rand_gen_seq_'+str(num)+'.fasta','fasta')
    
    print('BLAST running...........')
    ### BLAST Run
    subprocess.call([blast_bin_path+'/blastp','-query','./aug_data/rand_gen_seq/rand_gen_seq_'+str(num)+'.fasta',\
                 '-db','./aug_data/blast_db/blast_db',\
                 '-out','./aug_data/blast_out/rand_gen_seq_'+str(num)+'.xml',\
                 '-evalue',str(e_val),\
                 '-outfmt','5','-num_threads',str(blast_cpu_usage)]) 
    
    ### Read result xml file
    blast_qresult = SearchIO.parse('./aug_data/blast_out/rand_gen_seq_'+str(num)+'.xml','blast-xml')
    print('Check hit sequence.........')
    ### Check whether hit same sequence
    for i, q_res in enumerate(blast_qresult):
        switch_query = False
        for hsp in q_res.hsps:
            if hsp.hit.id in hsp.query.id:
                ## Check matched seq
                switch_query = True
                ## Matched seq count
                res = p.search(hsp.query.id)
                df_id.loc[df_id['ID']==res.group(1),'Count']+=1

        if switch_query:
             blast_matched_seq.append(rand_gen_seq_list[i])
    res = np.sum(df_id['Count']>=aug_pep_num)
    print(res, '/', df_id.shape[0])
    
    print('------------------------------')
SeqIO.write(blast_matched_seq,'./aug_data/augmented_peptides.fasta','fasta')