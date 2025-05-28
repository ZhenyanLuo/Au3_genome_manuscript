#https://github.com/google/deepconsensus/blob/r1.2/docs/quick_start.md install deepconsensus, (python=3.9.16)

#!/bin/bash
#PBS -q normal
#PBS -l mem=190GB
#PBS -l walltime=12:00:00
#PBS -l jobfs=2GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Apsidii_2023/
pbindex Au3.subreads.bam

#Save this one as .sh
n_total=500
subreads_path="/g/data/xf3/zl1602/Apsidii_2023/Au3.subreads.bam"
for ((n=1;n<=${n_total};n++));
do
function to_shard_id {
  echo "$( printf %03g "${1}")-of-$(printf "%03g" "${2}")"
}
shard_id="$(to_shard_id "${n}" "${n_total}")"
mkdir -p D01_${shard_id}
PBS_setting="#!/bin/bash
#PBS -q normal
#PBS -l mem=10GB
#PBS -l walltime=0:30:00
#PBS -l ncpus=16
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Apsidii_2023/Au3_${shard_id}
"
ccs_cmd="
ccs --min-rq=0.88 -j 16 --chunk=${n}/${n_total} ${subreads_path} Au3_${shard_id}.ccs.bam
pbindex Au3_${shard_id}.ccs.bam
"
echo -e "${PBS_setting}\n${ccs_cmd}" >Au3_${shard_id}.ccs.pbs.sh
done


#make actc scripts

n_total=500
subreads_path="/g/data/xf3/zl1602/Apsidii_2023/Au3.subreads.bam"
for ((n=1;n<=${n_total};n++));
do
function to_shard_id {
  echo "$( printf %03g "${1}")-of-$(printf "%03g" "${2}")"
}
shard_id="$(to_shard_id "${n}" "${n_total}")"
PBS_setting="#!/bin/bash
#PBS -q normal
#PBS -l mem=10GB
#PBS -l walltime=2:00:00
#PBS -l ncpus=3
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
cd /g/data/xf3/zl1602/Apsidii_2023/Au3_${shard_id}
"
actc_cmd="
actc ${subreads_path} Au3_${shard_id}.ccs.bam Au3_${shard_id}.subreads_to_ccs.bam
pbindex Au3_${shard_id}.subreads_to_ccs.bam
"
echo -e "${PBS_setting}\n${actc_cmd}" >Au3_${shard_id}/Au3_${shard_id}.actc.pbs.sh
done


#make deepconsensus script

n_total=500
for ((n=1;n<=${n_total};n++));
do
function to_shard_id {
  echo "$( printf %03g "${1}")-of-$(printf "%03g" "${2}")"
}
shard_id="$(to_shard_id "${n}" "${n_total}")"
PBS_setting="#!/bin/bash
#PBS -q normal
#PBS -l mem=60GB
#PBS -l walltime=6:00:00
#PBS -l ncpus=10
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3
#PBS -P xf3
set -xue
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate deepconsensus
"
id="Au3_${shard_id}"
subreads_path="/g/data/xf3/zl1602/Apsidii_2023/Au3.subreads.bam"
subreads_to_ccs="/g/data/xf3/zl1602/Apsidii_2023/${id}/${id}.subreads_to_ccs.bam"
ccs_bam="/g/data/xf3/zl1602/Apsidii_2023/${id}/${id}.ccs.bam"
model="/g/data/xf3/zl1602/Apsidii_2023/model"
PBS_d='$PBS_JOBFS'

deepconsensus_cmd="
cd ${PBS_d}
cp ${subreads_to_ccs} ${subreads_to_ccs}.pbi .
cp ${ccs_bam} ${ccs_bam}.pbi .
cp -r ${model} .
#Change batch_size and batch_zmws to half or less if running out of memory
deepconsensus run --subreads_to_ccs=${id}.subreads_to_ccs.bam --ccs_bam=${id}.ccs.bam --checkpoint=model/checkpoint --output=${id}.output.fastq --batch_size=2048 --batch_zmws=1000
cp *.output.* /g/data/xf3/zl1602/Apsidii_2023/Au3/${id}
"
echo -e "${PBS_setting}\n${deepconsensus_cmd}" >${id}/${id}.deepconsensus.pbs.sh
done


#Combine all output into one
find . -type f|grep 'output.fastq'|xargs mv output_fastq/
cd output_fastq
cat *.output.fastq >>deepconsensus_processed.fastq
