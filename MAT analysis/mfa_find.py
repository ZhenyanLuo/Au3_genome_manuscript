import sys
new_site = '/g/data/fa63/zl1602/software/SOG/lib/python3.12/site-packages'
sys.path.insert(0, new_site)
import pandas as pd
from Bio import SeqIO
from Bio import SearchIO
import os
import glob
from pathlib import Path
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from itertools import combinations
from pybedtools import BedTool

#Here import the fasta file
fna = list(SeqIO.parse('Au3_chr5.fna', 'fasta'))

#Find all the ORF  in the fragment sequences
orf_fwd = []

#Find all the ORF  in the fragment sequences
orf_fwd = []
for f in range(0, len(fna), 1):
    for i in range(0, len(fna[f].seq), 1):
        if fna[f].seq[i:i+3] == 'ATG':
            for j in range(i+3, len(fna[f].seq), 3):
                if fna[f].seq[j:j+3] in ['TAA', 'TAG', 'TGA']:
                    #if the ORF is too short, skip
                    if j-i < 120:
                        break
                    else:
                        header = fna[f].id + '_' + str(i) + '_' + str(j) + '_fwd'
                        orf_fwd.append(SeqRecord(fna[f].seq[i:j+3], id=header, description='orf'))
                    break


orf_rev = []
fna_rev = []
for f in range(0, len(fna), 1):
    fna_rev.append(SeqRecord(fna[f].reverse_complement().seq, id=fna[f].id, description='orf'))
for f in range(0, len(fna_rev), 1):
    for i in range(0, len(fna_rev[f].seq), 1):
        if fna_rev[f].seq[i:i+3] == 'ATG':
            for j in range(i+3, len(fna_rev[f].seq), 3):
                if fna_rev[f].seq[j:j+3] in ['TAA', 'TAG', 'TGA']:
                    #if the ORF is too short, skip
                    if j-i < 120:
                        break
                    else:
                        header = fna_rev[f].id + '_' + str(i) + '_' + str(j) + '_rev'
                        orf_rev.append(SeqRecord(fna_rev[f].seq[i:j+3], id=header, description='orf'))
                    break

orf_protein = []
for i in orf_fwd:
    orf_protein.append(SeqRecord(i.seq.translate(), id=i.id))
for i in orf_rev:
    orf_protein.append(SeqRecord(i.seq.translate(), id=i.id))

mfa_candidate_df = pd.DataFrame(columns=['aa', 'id', 'len'])
for i in orf_protein:
   pro_seq = i.seq
   #convert the pro_seq to string
   pro_seq = str(pro_seq)
   tmp = pd.DataFrame({'aa': pro_seq, 'id': i.id, 'len': len(pro_seq)}, index=[0])
   mfa_candidate_df = pd.concat([mfa_candidate_df, tmp])
mfa_candidate_df.to_csv('mfa_candidate.csv', index=False, header=True, sep='\t')
aliphatic_residues = ['A', 'V', 'I', 'L', 'M']
mfa_candidate_df_2 = mfa_candidate_df[mfa_candidate_df['aa'].str.contains(r'C[AVLI]{2}[A-Z]\*')]
mfa_candidate_df_3 = mfa_candidate_df[mfa_candidate_df['aa'].str.contains('SHIC')]
mfa_candidate_df_2.to_csv('mfa_candidate_2.csv', index=False, header=True, sep='\t')
mfa_candidate_df_3.to_csv('mfa_candidate_3.csv', index=False, header=True, sep='\t')