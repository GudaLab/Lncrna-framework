#!/bin/bash
#SBATCH --job-name=star # Job name
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
module load star/2.7.11b
ulimit -n 80000
call="STAR --runThreadN 50 --genomeDir human_genome/STAR_index --readFilesIn ${sample}_1_clean.fastq.gz ${sample}_2_clean.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}_ --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif"
echo $call
eval $call
end=`date +%s`
runtime=$((end-start))
echo $runtime