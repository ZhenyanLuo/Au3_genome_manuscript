import sys
new_site = '/g/data/fa63/zl1602/software/SOG/lib/python3.12/site-packages'
sys.path.insert(0, new_site)
import pandas as pd
from pybedtools import BedTool
import glob
import os

#Assign a number for running the script in the cluster
CpG_sites = pd.read_csv("merged_CpG_site.bed", sep="\t", header=None)
CpG_sum = CpG_sites.groupby([0]).size().reset_index(name='counts')
CpG_sum.columns = ['chrom', 'CpG_counts']
CpG_sum_2 =sum(CpG_sum['CpG_counts'])
TE_bed = pd.read_csv("Au3_TE_merged.bed", sep="\t", header=None)
#Convert the TE bed file to bedtools object
TE_bed = BedTool.from_dataframe(TE_bed)
#Sort the bed file
TE_bed = TE_bed.sort()

CDS_bed = pd.read_csv("Au3_CDS.bed", sep="\t", header=None)
CDS_bed = BedTool.from_dataframe(CDS_bed)
CDS_bed = CDS_bed.sort()

UTR3_bed = pd.read_csv("Au3_UTR3.bed", sep="\t", header=None)
#Only keep rows with 'T1'
UTR3_bed = UTR3_bed[UTR3_bed[3].str.contains('T1')]
#Only keep the frist 3 columns
UTR3_bed = UTR3_bed.iloc[:, 0:3]
UTR3_bed = BedTool.from_dataframe(UTR3_bed)
UTR3_bed = UTR3_bed.sort()

UTR5_bed = pd.read_csv("Au3_UTR5.bed", sep="\t", header=None)
#Only keep rows with 'T1'
UTR5_bed = UTR5_bed[UTR5_bed[3].str.contains('T1')]
#Only keep the frist 3 columns
UTR5_bed = UTR5_bed.iloc[:, 0:3]
UTR5_bed = BedTool.from_dataframe(UTR5_bed)
UTR5_bed = UTR5_bed.sort()


#assign the label to the output file, based on the value passed in the command line
run_label = sys.argv[1]
#Create a random seed for the shuffle command
import random
#assign the random seed to the shuffle command
random_seed = random.randint(1, 1000000)
cmd = f"bedtools shuffle -i Au3_5mCpG_filtered_3col.bed -g Au3.genome.size.fixed  -chrom -noOverlapping -incl merged_CpG_site.bed -seed {random_seed} > shuffle_{run_label}.bed"
os.system(cmd)
shuffled_bed = pd.read_csv(f"shuffle_{run_label}.bed", sep="\t", header=None)
shuffled_bed = shuffled_bed.dropna()
shuffled_bed = BedTool.from_dataframe(shuffled_bed)
shuffled_bed = shuffled_bed.sort()
intersect_TE = TE_bed.coverage(shuffled_bed)
intersect_TE_bed = intersect_TE.to_dataframe()
cov_sum = sum(intersect_TE_bed['name'])
TE_simulated_ratio = cov_sum/CpG_sum_2
intersect_CDS = CDS_bed.coverage(shuffled_bed)
intersect_CDS_bed = intersect_CDS.to_dataframe()
cov_sum = sum(intersect_CDS_bed['name'])
CDS_simulated_ratio = cov_sum/CpG_sum_2
intersect_UTR3 = UTR3_bed.coverage(shuffled_bed)
intersect_UTR3_bed = intersect_UTR3.to_dataframe()
cov_sum = sum(intersect_UTR3_bed['name'])
UTR3_simulated_ratio = cov_sum/CpG_sum_2
intersect_UTR5 = UTR5_bed.coverage(shuffled_bed)
intersect_UTR5_bed = intersect_UTR5.to_dataframe()
cov_sum = sum(intersect_UTR5_bed['name'])
UTR5_simulated_ratio = cov_sum/CpG_sum_2
final_result = pd.DataFrame({'run_label': [run_label], 'TE_simulated_ratio': [TE_simulated_ratio], 'CDS_simulated_ratio': [CDS_simulated_ratio], 'UTR3_simulated_ratio': [UTR3_simulated_ratio], 'UTR5_simulated_ratio': [UTR5_simulated_ratio]})
final_result.to_csv(f"shuffle_result2/{run_label}.txt", sep="\t", index=False, header=False)
#delete the shuffled file
os.remove(f"shuffle_{run_label}.bed")