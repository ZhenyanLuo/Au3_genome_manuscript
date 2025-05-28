import pandas as pd
import multiprocessing as mp
from functools import partial
from pathlib import Path
import re
import pandas as pd
from pathlib import Path
import numpy as np
from Bio import SeqIO

fna_hapA = list(SeqIO.parse('../../Au3_CHR_hapA_v3.fasta', 'fasta'))
fna_hapB = list(SeqIO.parse('../../Au3_CHR_hapB_v3.fasta', 'fasta'))

fna_rename = []
for i in range(0, len(fna_hapA)):
    seq = fna_hapA[i].seq
    if re.search(r'HapA', fna_hapA[i].id):
        hap = 'A'
    else:
        hap = 'B'
    chr = fna_hapA[i].id.split('_')[3]
    chr = chr.replace('CHR', 'Chr')
    header = chr+hap
    if fna_hapA[i].id == 'APSI_AU3_HapA_CHR14_ab':
        header = 'Chr14B'
    fna_rename.append(SeqIO.SeqRecord(seq, id=header, description=''))
for i in range(0, len(fna_hapB)):
    seq = fna_hapB[i].seq
    if re.search(r'HapA', fna_hapB[i].id):
        hap = 'A'
    else:
        hap = 'B'
    chr = fna_hapB[i].id.split('_')[3]
    chr = chr.replace('CHR', 'Chr')
    header = chr+hap
    fna_rename.append(SeqIO.SeqRecord(seq, id=header, description=''))

bed_5mC = pd.read_table('../Au3_5mCpG.bed', header=None)
bed_5mC.columns = ['chr', 'start', 'end', 'score','reads_count']

def process_item(item):
    chr = item.id
    tmp = bed_5mC[bed_5mC['chr'] == chr]
    seq = item.seq
    result = pd.DataFrame()
    result0 = pd.DataFrame()
    for index, row in tmp.iterrows():
        start = row['start']
        end = row['end'] + 1
        if seq[start:end] == 'CG':
            result = pd.concat([result, pd.DataFrame([row])])
        else:
            result0 = pd.concat([result0, pd.DataFrame([row])])
    file_name = f'Au3_5mCpG_{chr}.bed'
    file_name0 = f'Au3_5mC_{chr}.bed'
    result.to_csv(file_name, sep='\t', header=False, index=False)
    result0.to_csv(file_name0, sep='\t', header=False, index=False)
    return

def process_files_parallel(seqname, num_threads):
    pool = mp.Pool(num_threads)
    pool.map(process_item, seqname)
    pool.close()
    pool.join()

if __name__ == '__main__':
    process_files_parallel(fna_rename, 18)

