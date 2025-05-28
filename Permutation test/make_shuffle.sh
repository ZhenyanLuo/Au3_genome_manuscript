for i in {1..1000} :
do
    cat << EOF > "shuffle_${i}.pbs"
#!/bin/bash
#PBS -q normal
#PBS -l mem=50GB
#PBS -l walltime=00:20:00
#PBS -l jobfs=180GB
#PBS -l ncpus=1
#PBS -P xf3
#PBS -l wd
#PBS -l storage=scratch/fa63+gdata/fa63+gdata/xf3+scratch/xf3
module load python3
module load bedtools
cd /scratch/xf3/zl1602/permutate
python3 20250213_shuffle.py ${i}
EOF
    echo "Created PBS script: shuffle_${i}.pbs"
done

