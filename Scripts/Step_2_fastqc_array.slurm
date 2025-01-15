#!/bin/bash
#SBATCH --job-name=fastqc # Job name
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
module load fastqc/0.12.1
call="fastqc ${sample}_1.fastq.gz ${sample}_2.fastq.gz -o qc_results"
echo $call
eval $call
end=`date +%s`
runtime=$((end-start))
echo $runtime