#Au3
funannotate train -i APSI_v3.hapA.fasta.masked -l APSI_v3_A_R1.fastq.gz -r APSI_v3_A_R2.fastq.gz --memory 200G --jaccard_clip --cpus 24 --species Austropuccinia_psidii --strain Au3 -o Au3_hapA_train \
  --nanopore_mrna APSI_v3_A_ONT_mapped.fastq.gz
funannotate predict -i APSI_v3.hapA.fasta.masked -o Au3_hapA_train -s "Austropuccinia_psidii" --strain Au3_hapA --cpus 12 --force --name MK676
funannotate update -i Au3_hapA_train --cpus 24 --pasa_db mysql --jaccard_clip

funannotate train -i APSI_v3.hapB.fasta.masked -l APSI_v3_B_R1.fastq.gz -r APSI_v3_B_R2.fastq.gz --memory 200G --jaccard_clip --cpus 24 --species Austropuccinia_psidii --strain Au3 -o Au3_hapB_train \
  --nanopore_mrna APSI_v3_B_ONT_mapped.fastq.gz
funannotate predict -i APSI_v3.hapB.fasta.masked -o Au3_hapB_train -s "Austropuccinia_psidii" --strain Au3_hapB --cpus 12 --protein_evidence hapA_Au3.faa $FUNANNOTATE_DB/uniprot_sprot.fasta --force --name MK675
funannotate update -i Au3_hapB_train --cpus 24 --pasa_db mysql --jaccard_clip
#!/bin/bash
#PBS -q normal
#PBS -l mem=80GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=16
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/fa63+gdata/if89+gdata/fa63
#PBS -P fa63
set -xue
module load java/jdk-13.33
module load python3
export PATH=$PATH:/g/data/fa63/share/interproscan-5.61-93.0
cd /g/data/fa63/zl1602/Annotation
interproscan.sh -i Austropuccinia_psidii_Au3_hapA.proteins.fa -dp -cpu 16 --output-dir Au3_hapA_ipr

#!/bin/bash
#PBS -q normal
#PBS -l mem=80GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=16
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -l storage=scratch/fa63+gdata/if89+gdata/fa63
#PBS -P fa63
set -xue
module load java/jdk-13.33
module load python3
export PATH=$PATH:/g/data/fa63/share/interproscan-5.61-93.0
cd /g/data/fa63/zl1602/Annotation
interproscan.sh -i Austropuccinia_psidii_Au3_hapB.proteins.fa -dp -cpu 16 --output-dir Au3_hapB_ipr


for i in `cat list`;
do
setting="#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=4GB
#SBATCH --job-name=${i}
#SBATCH --requeue
source /opt/conda/etc/profile.d/conda.sh
conda activate /mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate
export PASAHOME=/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/opt/pasa-2.5.3
export TRINITY_HOME=/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/opt/trinity-2.8.5
export EVM_HOME=/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/opt/evidencemodeler-1.1.1
export AUGUSTUS_CONFIG_PATH=/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/config/
export FUNANNOTATE_DB=/mnt/data/dayhoff/home/scratch/u6575017/DB
export PERL5LIB=/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/lib/perl5/site_perl/:/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/lib/perl5/core_perl/:/mnt/data/dayhoff/home/u6575017/.conda/envs/funannotate/lib/perl5/vendor_perl/:/home/u6575017/.conda/envs/funannotate/lib/perl5/5.32/site_perl:/mnt/data/dayhoff/home/u6575017/software/signalp-4.1/lib/
export PATH=/usr/bin/:/home/u6575017/.conda/envs/funannotate/bin/:/mnt/data/dayhoff/home/u6575017/software/eggnog-mapper_2.1.9/eggnog-mapper-2.1.9/:/mnt/data/dayhoff/home/u6575017/software/signalp-4.1
export GENEMARK_PATH=/mnt/data/dayhoff/home/u6575017/software/gmes_linux_64_4
set -xue
"
prefix=`echo ${i}|sed -e 's/_hap.*//'`
cmd="
cd /mnt/data/dayhoff/home/scratch/u6575017/Annotation/${i}
funannotate annotate -i ${i}_train/update_results -o ${i}_annotate --cpus 15 --iprscan ${i}_interpro.xml --strain ${i} --busco_db basidiomycota --sbt ${prefix}_template.sbt
"
echo -e "${setting}\n${cmd}" >${i}_annotate.pbs.sh
done



