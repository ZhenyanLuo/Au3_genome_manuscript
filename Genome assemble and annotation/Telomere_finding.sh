#!/bin/bash
#PBS -q normal
#PBS -l mem=50GB
#PBS -l walltime=3:00:00
#PBS -l ncpus=1ls
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -l storage=scratch/xf3+gdata/xf3+gdata/fa63+scratch/fa63
#PBS -P fa63
export PATH=$PATH:$HOME/.cargo/bin
cd /g/data/fa63/zl1602/SOG/genome
for i in $(ls *.fa); do
    prefix=`echo $i | cut -d '.' -f 1`
    tidk search --string 'AACCCT' --output ${i}_tidk --dir ${i}_tidk -w 50 ${i}
done

for i in *_tidk ; do file=`ls ${i}` ; mv ${i}/${file} tidk/ ; done




for i in *_tidk.tsv ;
do
prefix=`echo ${i}|cut -d'.' -f1`
tidk plot --tsv ${i} 
mv tidk-plot.svg ${prefix}.svg
done



#Only keep the lines which have more than 4 match in column 3 or column 4
for i in *_tidk.tsv; 
do
    prefix=`echo $i | cut -d '.' -f 1 |sed -e 's/Austropuccinia_psidii/Apsi/' `
    cat $i | awk -v OFS='\t' '{if($3>=4 || $4>=4) print $0}' > ${prefix}_f.tsv
done
