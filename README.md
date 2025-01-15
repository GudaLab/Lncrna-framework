# lncRNA Analysis Pipeline
LncRNA framework for identification of lncRNA dysregulations.

## Overview
LncRNA framework pipeline facilitates comprehensive lncRNA analysis, including preprocessing, alignment, transcript quantification, differential expression analysis, and functional annotation. The tools and steps outlined in this document ensure accurate and reproducible results. The below steps are for processing a single sample. Please see the scripts folder for batch processing.

---

## Prerequisites

### Required Tools
1. [FastQC](https://github.com/s-andrews/FastQC)
2. [fastp](https://github.com/OpenGene/fastp)
3. [STAR](https://github.com/alexdobin/STAR)
   - Use `--quantMode GeneCounts` for quantification.
4. [SAMtools](http://www.htslib.org/)
5. [StringTie](https://github.com/gpertea/stringtie)
   - Use `-e -B` for quantification.
   - Use `StringTie Merge` with `-G` for reference-based merging.
6. [gffcompare](https://github.com/gpertea/gffcompare) (use `-G` for compatibility).
7. [awk](https://www.gnu.org/software/gawk/) (for processing files).
8. [CPAT](https://pypi.org/project/CPAT/) (for coding potential assessment).
9. [featureCounts](https://subread.sourceforge.net/featureCounts.html) (use `-p -B -C` for paired-end data).
10. [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) (for differential expression analysis).
11. [FEELnc](https://github.com/tderrien/FEELnc) (for functional annotation).
12. [Miranda](https://github.com/hacktrackgnulinux/miranda) (for target prediction).
13. [ClusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) (for pathway enrichment).
14. Visualization tools:
    - [ggplot2](https://ggplot2.tidyverse.org/)
    - [igraph](https://igraph.org/)
15. [GENCODE Annotations](https://www.gencodegenes.org/human/)

### Required Reference
- **Genome:** GRCh38.p14
- **Annotation Files:** GENCODE v47 GTF/GFF

---

## Pipeline Steps

### Step 1: Download Data from SRA
```bash
fastq-dump SRR24709142 --split-3 --gzip
```

### Step 2: Quality Control of Raw Reads
```bash
fastqc SRR24709142_1.fastq.gz SRR24709142_2.fastq.gz -o qc_results
```

### Step 3: Trimming and Filtering
```bash
fastp -i SRR24709142_1.fastq.gz -I SRR24709142_2.fastq.gz -o SRR24709142_1_clean.fastq.gz -O SRR24709142_2_clean.fastq.gz -h SRR24709142_report.html --trim_poly_g --cut_mean_quality 20 --length_required 50
```

### Step 4: Generate STAR Index
```bash
STAR --runThreadN 50 --runMode genomeGenerate --genomeDir human_genome/STAR_index --genomeFastaFiles human_genome/GRCh38.p14.genome.fa --sjdbGTFfile human_genome/gencode.v47.long_noncoding_RNAs.gtf
```

### Step 5: Align Reads to Reference Genome
```bash
STAR --runThreadN 50 --genomeDir human_genome/STAR_index --readFilesIn SRR24709142_1_clean.fastq.gz SRR24709142_2_clean.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix SRR24709142_ --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif
```

### Step 6: Assemble Transcripts
```bash
cufflinks -g human_genome/gencode.v47.long_noncoding_RNAs.gtf -o SRR24709142_cufflinks_out -p 8 SRR24709142_Aligned.sortedByCoord.out.bam
```

### Step 7: Create Transcript Merge List
Create a `mergelist.txt` file with paths to all `transcripts.gtf` files.

### Step 8: Merge Transcripts
```bash
awk '!seen[$1, $4, $5]++' merged.gtf > merged_without_duplicates.gtf
```

### Step 9: Remove Mitochondrial Chromosome
```bash
awk '/^>/{p=($0!="chrM MT")} p' human_genome/GRCh38.p14.genome.fa > human_genome/GRCh38.p14.genome_no_chrM.fa
```

### Step 10: Assess Coding Potential
Download pre-built models and run:
```bash
cpat.py -g filtered_lncRNAs.fa -x Human_Hexamer.tsv -d Human_logitModel.RData -o cpat_output.txt
```

### Step 11: Quantify Transcript Abundance
```bash
featureCounts -a merged_without_duplicates.gtf -o SRR24709142_counts.txt -T 8 -p -B -C SRR24709142_Aligned.sortedByCoord.out.bam
```

### Step 12: Process Count Files
Extract and combine processed data:
```bash
for file in *_counts.txt; do
  sed '1d' "$file" | awk -v OFS="\t" '{print $1, $7}' > "processed_$file"
done
cut -f1 processed_SRR24709142_counts.txt > combined_file.txt
for file in processed_*_counts.txt; do
  awk '{print $2}' "$file" | paste combined_file.txt - > temp.txt
  mv temp.txt combined_file.txt
done
```

### Step 13: Differential Expression Analysis
Use R with DESeq2 for differential expression analysis and visualization.

### Step 14: Functional Annotation
```bash
FEELnc_filter.pl -i filtered_gencode.gtf -a gencode.v47.annotation.gtf -b transcript_biotype=protein_coding -l > candidate_lncRNA.gtf
awk 'BEGIN{FS="\t"; OFS="\t"} {match($9, /gene_id \"([^\"]+)\"/, a); if (a[1] != "") print a[1]}' candidate_lncRNA.gtf | sort | uniq > lncRNA_target_genes.txt
```

### Step 15: Pathway Enrichment Analysis
Use `clusterProfiler` in R for pathway enrichment analysis.

---

## Output Files
1. **Processed Counts**: `combined_file.txt`
2. **Differential Expression Results**: `differential_expression_results_with_labels.csv`
3. **Pathway Enrichment Results**: `pathway_enrichment_results_ENSEMBL_ID.csv`

---

## Notes
- Ensure all tools are installed and configured in a dedicated Conda environment.
- Use sufficient computational resources for STAR alignment and StringTie assembly.
- Verify input file formats before running each step.
