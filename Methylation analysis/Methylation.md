#Install  pbmm2-1.13.1 (onto xf3/mafft env)
conda install bioconda::pbmm2 -y

#Index reference 
pbmm2 index Au3.fasta Au3.mmi --preset CCS



#Align ccs to reference, preset CCS was used
Alignment modes of --preset:
    SUBREAD     : -k 19 -w 19    -o 5 -O 56 -e 4 -E 1 -A 2 -B 5 -z 400 -Z 50  -r 2000   -g 5000
    CCS or HiFi : -k 19 -w 19 -u -o 6 -O 26 -e 2 -E 1 -A 1 -B 4 -z 400 -Z 50  -r 2000   -g 5000
    ISOSEQ      : -k 15 -w 5  -u -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 200000 -g 2000 -C 5 -G 200000
    UNROLLED    : -k 15 -w 15    -o 2 -O 32 -e 1 -E 0 -A 1 -B 2 -z 200 -Z 100 -r 2000   -g 10000

#!/bin/bash
#PBS -q hugemem
#PBS -l mem=800GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
module load java/jdk-8.40
module load python2/2.7.17
module load parallel/20191022
module load bwa/0.7.17
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
cd /scratch/xf3/zl1602/Methylation
conda activate mafft
pbmm2 align Au3.fasta m64123_201211_162803.Q20.fastq Au3.movie.bam --sort --rg '@RG\tID:myid\tSM:mysample' --preset CCS --num-threads 32



#Install pbCpG
wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
tar -xzf pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz


#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
cd /scratch/xf3/zl1602/Methylation
pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam Au3.movie.bam \
  --output-prefix Au3.movie.pbmm2 \
  --model pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 32







###Also try ccsmeth
#Install env
conda env create --prefix /scratch/xf3/zl1602/ccsmeth_env -f ccsmeth/environment.yml

conda activate /scratch/xf3/zl1602/ccsmeth_env




#-------------------------------------Ignore Above Command

#!/bin/bash
#PBS -q hugemem
#PBS -l mem=300GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=600GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation
ccs m64123_201211_162803.subreads.bam Au3.hifi_reads.bam --hifi-kinetics -j 32


#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=1
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation
pbindex m64123_201211_162803.subreads.bam


n_total=500
subreads_path="/g/data/xf3/zl1602/Methylation/m64123_201211_162803.subreads.bam"
for ((n=1;n<=${n_total};n++));
do
function to_shard_id {
  echo "$( printf %03g "${1}")-of-$(printf "%03g" "${2}")"
}
shard_id="$(to_shard_id "${n}" "${n_total}")"
PBS_setting="#!/bin/bash
#PBS -q normal
#PBS -l mem=15GB
#PBS -l walltime=1:30:00
#PBS -l ncpus=16
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
#PBS -l jobfs=200GB
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
module load parallel/20191022
cd /g/data/xf3/zl1602/Methylation/Au3_chunk
"
ccs_cmd="
ccs m64123_201211_162803.subreads.bam Au3.ccs.${shard_id}.bam --chunk ${n}/${n_total} --hifi-kinetics -j 16
"
echo -e "${PBS_setting}\n${ccs_cmd}" >Au3_${shard_id}.ccs.pbs.sh
done



#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation/Au3_chunk
module load samtools
samtools merge -@ 32 Au3.ccs.bam Au3.ccs.*-of-500.





#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation/Au3_chunk
jasmine Au3.ccs.bam Au3.hifi_reads.bam -j 32


###---For no chunked data (the ccs step takes longer time but the file size is bigger than chunked one so let's keep both files)
#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation/
jasmine Au3.hifi_reads.bam Au3.hifi.movie.bam -j 32
###----





#------------create index
pbmm2 index Au3.fasta Au3.mmi --preset CCS
#------------


#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation/Au3_chunk
pbmm2 align Au3.mmi Au3.hifi_reads.bam Au3.hifi_pbhmm2.bam --preset CCS --sort -j 32 

#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation/Au3_chunk
pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam Au3.hifi_pbhmm2.bam \
  --output-prefix Au3.hifi_chunked.pbhmm2 \
  --model pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 32





########
###---For no chunked data (the ccs step takes longer time but the file size is bigger than chunked one so let's keep both files)
#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation/
pbmm2 index Au3.fasta Au3.mmi --preset CCS
pbmm2 align Au3.mmi Au3.hifi.movie.bam Au3.hifi.movie_pbhmm2.bam --preset CCS --sort -j 32 


#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Methylation
pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam Au3.hifi.movie_pbhmm2.bam \
  --output-prefix Au3.hifi.movie_nochunked.pbhmm2 \
  --model pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 32



cat Au3.hifi.movie_nochunked.pbhmm2.combined.bed |awk -v OFS='\t' '{if($4 >=70 )print $0}' >Methylation_Analysis/Au3.filtered.CpG.bed
sed -i 's/APSI_//' Au3.filtered.CpG.bed
cat Au3.filtered.CpG.bed|awk -v OFS='\t' '{print $1,$2,$3}'|sort -k1,1 -k2,2n >Au3.CpG.bed
samtools faidx Au3.fasta.fai
cat Au3.fasta.fai|awk -v OFS='\t' '{print $1, $2}'|sort -k1,1 >Au3.fasta.genome
bedtools makewindows -w 10000 -g Au3.fasta.genome > Au3.fasta.window
bedtools coverage -a Au3.fasta.window -b Au3.CpG.bed -counts >Au3.CpG.cov.bed





#--------------Now make for MF1
#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=16
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P fa63
set -xue
cd /g/data/xf3/zl1602/Methylation/MF1
tar -I 'pigz -p 16' -xvf XBSAU.20221222_134155.PACBIO_DATA.2-of-3.tgz



#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=2
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P fa63
set -xue
cd /g/data/xf3/zl1602/Methylation/MF1/MF1/PACBIO_DATA
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/fa63/share/miniconda3/envs/Pacbio
pbindex MF1.subreads.bam




#!/bin/bash
#PBS -q hugemem
#PBS -l mem=300GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=600GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P fa63
set -xue
cd /g/data/xf3/zl1602/Methylation/MF1/MF1/PACBIO_DATA
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/fa63/share/miniconda3/envs/Pacbio
ccs MF1.subreads.bam MF1.hifi_reads.bam --hifi-kinetics -j 32




#!/bin/bash
#PBS -q normal
#PBS -l mem=120GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=200GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/fa63/share/miniconda3/envs/Pacbio
cd /g/data/xf3/zl1602/Methylation/MF1/MF1/PACBIO_DATA
jasmine MF1.hifi_reads.bam MF1.hifi.movie.bam -j 32

#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/fa63/share/miniconda3/envs/Pacbio
cd /g/data/xf3/zl1602/Methylation/MF1/MF1/PACBIO_DATA
pbmm2 index MF1.fasta MF1.mmi --preset CCS
pbmm2 align MF1.mmi MF1.hifi.movie.bam MF1.hifi.movie_pbhmm2.bam --preset CCS --sort -j 32 



#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=32
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/fa63/share/miniconda3/envs/Pacbio
cd /g/data/xf3/zl1602/Methylation/MF1/MF1/PACBIO_DATA
/g/data/xf3/zl1602/Methylation/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam MF1.hifi.movie_pbhmm2.bam \
  --output-prefix MF1.hifi.movie.pbhmm2 \
  --model /g/data/xf3/zl1602/Methylation/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 32





#!/bin/bash
#PBS -q normal
#PBS -l mem=150GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=2
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/Au3_ONT
module load python3
source /g/data/fa63/zl1602/software/SOG/bin/activate
export PATH=$PATH:/home/106/zl1602/.local/bin/
pod5 convert fast5 fast5_pass/*.fast5 --output PAO37085_Au3.pod5


#!/bin/bash
#PBS -q normal
#PBS -l mem=150GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=2
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/Au3_ONT
module load python3
source /g/data/fa63/zl1602/software/SOG/bin/activate
export PATH=$PATH:/home/106/zl1602/.local/bin/
pod5 convert fast5 /g/data/fa63/bxs800/ONT_gDNA87_469_BS/20230502_Au_3/20230502_1711_2D_PAO37085_0d6b041c/fast5_fail/*.fast5 --output /g/data/fa63/zl1602/PAO37085_Au3_fail.pod5






pod5 inspect debug

dorado download --model dna_r10.4.1_e8.2_260bps_sup@v4.1.0 

export LD_LIBRARY_PATH=<PATH_TO_DORADO>/dorado-x.y.z-linux-x64/lib:$LD_LIBRARY_PATH


#!/bin/bash
#PBS -q gpuvolta
#PBS -l mem=200GB
#PBS -l walltime=24:00:00
#PBS -l ngpus=4
#PBS -l ncpus=48
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/Au3_ONT
module load cuda/11.4.1
module load openmpi/4.1.0
module load python3
module load minimap2
source /g/data/fa63/zl1602/software/SOG/bin/activate
export PATH=$PATH:/home/106/zl1602/.local/bin/:/g/data/xf3/miniconda/share_software/dorado-0.7.3-linux-x64/bin/
export LD_LIBRARY_PATH=/g/data/xf3/miniconda/share_software/dorado-0.7.3-linux-x64/lib:$LD_LIBRARY_PATH
dorado basecaller dna_r10.4.1_e8.2_260bps_sup@v4.1.0  PAO37085_Au3.pod5   --modified-bases 5mCG_5hmCG  --device "cuda:all" > /scratch/xf3/zl1602/PAO37085_Au3.basecalled.bam

#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=24
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /scratch/xf3/zl1602/Au3_ONT_Methy
module load python3
module load minimap2
source /g/data/fa63/zl1602/software/SOG/bin/activate
export PATH=$PATH:/home/106/zl1602/.local/bin/:/g/data/xf3/miniconda/share_software/dorado-0.7.3-linux-x64/bin/
export LD_LIBRARY_PATH=/g/data/xf3/miniconda/share_software/dorado-0.7.3-linux-x64/lib:$LD_LIBRARY_PATH
dorado aligner -t 24 Au3.mmi PAO37085_Au3.basecalled.bam >PAO37085_Au3.alinged.bam



#!/bin/bash
#PBS -q normal
#PBS -l mem=80GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=24
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /scratch/xf3/zl1602/Au3_ONT_Methy
module load samtools
samtools sort  PAO37085_Au3.alinged.bam -@ 24 -o  PAO37085_Au3.sorted.bam
samtools index -@ 24 PAO37085_Au3.sorted.bam



#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=24
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /scratch/xf3/zl1602/Au3_ONT_Methy
module load python3
module load minimap2
export PATH=$PATH:/home/106/zl1602/.local/bin/:/g/data/xf3/miniconda/share_software/dorado-0.7.3-linux-x64/bin/:/g/data/xf3/miniconda/share_software/modkit
modkit pileup  -t 24 -r Au3.fna  PAO37085_Au3.sorted.bam PAO37085_Au3.methylation_calls.bed 
modkit extract   -t 24 -r Au3.fna  PAO37085_Au3.sorted.bam PAO37085_Au3_methylation_per_read.tsv
modkit pileup -t 24 -r Au3.fna  PAO37085_Au3.sorted.bam --bedgraph modkit/


#---------------------------Now look at 6mA
#We already have the kinetics hifi reads 

#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=24
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/6mA
export PATH=$PATH:/home/106/zl1602/.cargo/bin/
ft predict-m6a -t 24 Au3.hifi_reads.bam Au3.hifi_reads_m6a.bam

#!/bin/bash
#PBS -q hugemem
#PBS -l mem=200GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=24
#PBS -l jobfs=200GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/6mA
export PATH=$PATH:/home/106/zl1602/.cargo/bin/
module load minimap2
module load samtools
minimap2 -x map-hifi Au3.fna Au3.hifi_reads_m6a.bam -a -t 24 >Au3.hifi_reads_m6a.map.sam
samtools view -bS -@ 24 Au3.hifi_reads_m6a.map.sam|samtools sort -@ 24 -o Au3.hifi_reads_m6a.sort.bam -
samtools index -@ 24 Au3.hifi_reads_m6a.sort.bam




#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=6
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/6mA
module load python3
module load minimap2
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate Pacbio
pbindex Au3.hifi_reads.bam
#For add ip and wd
ccs-kinetics-bystrandify Au3.hifi_reads.bam Au3.hifi_reads_kinetics.bam -j 0

#!/bin/bash
#PBS -q hugemem
#PBS -l mem=300GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=32
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/6mA/actc
module load python3
module load minimap2
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate Pacbio
actc -j 32 m64123_201211_162803.subreads.bam Au3.hifi_reads.bam Au3.hifi.actc.bam
pbindex Au3.hifi.actc.bam


#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=48
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/zl1602/6mA/actc/test2
module load python3
module load minimap2
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate Pacbio
pbmm2 align Au3.mmi Au3.hifi_reads_kinetics.bam Au3.hifi_kinetics_pbmm2.bam --preset CCS --sort -J 24 -j 24




#!/bin/bash
#PBS -q hugemem
#PBS -l mem=200GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=32
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/miniconda/share_software/kineticsTools/kineticsTools
module load python3
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate py3.10.0
export PATH=/g/data/xf3/miniconda/share_software/kineticsTools/kineticsTools/:$PATH
wkdir=/g/data/xf3/zl1602/6mA/actc
python3 ipdSummary.py ${wkdir}/Au3.hifi.actc.bam --reference ${wkdir}/Au3.fasta --identify m6A,m4C,m5C_TET --numWorkers 32 --gff ${wkdir}/Au3.ipdSummary.gff


wkdir=/g/data/xf3/zl1602/6mA/actc/test2/
python3 ipdSummary.py ${wkdir}/Au3.hifi_kinetics_pbmm2.bam --reference ${wkdir}/Au3.fasta --identify m6A,m4C,m5C_TET --numWorkers 48 --gff ${wkdir}/Au3.ipdSummary.gff


#!/bin/bash
#PBS -q normal
#PBS -l mem=120GB
#PBS -l walltime=24:00:00
#PBS -l ncpus=48
#PBS -l jobfs=300GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate Pacbio
export PATH=$PATH:/home/106/zl1602/.cargo/bin/
module load pytorch/1.10.0
cd /g/data/xf3/zl1602/6mA/actc/ft
pbmm2 align Au3.fasta Au3.hifi_reads.bam Au3.hifi_pbmm2.bam --preset CCS --sort -J 0 -j 0
ft predict-m6a -t 48 Au3.hifi_pbmm2.bam Au3.hifi_kinetics_pbmm2.m6A.bam




#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=48
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
module load python3
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate Pacbio
export PATH=/g/data/xf3/miniconda/share_software/kineticsTools/kineticsTools/:$PATH
wkdir=/g/data/xf3/zl1602/6mA/actc
cd ${wkdir}/test_pbmm2
pbindex Au3.hifi.actc.bam
pbmm2 align --preset CCS Au3.fasta Au3.hifi.actc.bam  Au3.hifi.actc_to_ref.bam -j 24 -J 24  --sort 



#Get 6mA and 4mC result

for i in `cat list`;
do
chr=${i}
cmd="#!/bin/bash
#PBS -q normal
#PBS -l mem=180GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=24
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /g/data/xf3/miniconda/share_software/kineticsTools/kineticsTools
module load python3
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate py3.10.0
wkdir=/g/data/xf3/zl1602/6mA/actc/test2
python3 ipdSummary.py ${wkdir}/Au3.hifi_kinetics_pbmm2.bam --reference ${wkdir}/Au3.fasta --identify m6A,m4C --numWorkers 24 --gff ${wkdir}/chunk/${chr}.ipdSummary.gff --csv ${wkdir}/chunk/${chr}.ipdSummary.csv --methylFraction --refContigs ${i}
"
echo -e "$cmd" >${i}.ipdsummary.pbs.sh
done
















import sys
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

candidate=list(SeqIO.parse('mfa_WG.fna', 'fasta'))


candidate_orf_aa =[]

for record in candidate:
	for i in range(0, len(record.seq), 1):
		if record.seq[i:i+3] == 'ATG':
			for j in range(i+3, len(record.seq), 3):
				if record.seq[j:j+3] in ['TAA', 'TAG', 'TGA']:
					if len(record.seq[i:j+3].translate()) >=30:
						candidate_orf_aa.append(SeqRecord(record.seq[i:j+3].translate(), id=record.id))