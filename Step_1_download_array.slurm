#!/bin/bash
#SBATCH --job-name=ScRNA # Job name
#SBATCH --ntasks=10
#SBATCH --time=20-23:00:00
#SBATCH --mem=50Gb # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=1-57
#SBATCH --output=slurmout/%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/%A_%a.err # File to which STDERR will be written
start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" filename.txt`
echo "SAMPLE: ${sample}"
module load sra-tools/3.1.1
call="fastq-dump --split-3 --gzip ${sample}"
echo $call
eval $call
end=`date +%s`
runtime=$((end-start))
echo $runtime