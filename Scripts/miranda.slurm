#!/bin/bash
#SBATCH --job-name=miranda # Job name
#SBATCH --ntasks=25
#SBATCH --time=20-23:00:00
#SBATCH --mem=200Gb # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=miranda1.out
#SBATCH --error=miranda1.err
module load miranda/3.3a
miranda filtered_lncRNAs_differential_genes.fa gencode.v47.pc_transcripts.fa -sc 140 -en -20 -out miranda_output.txt