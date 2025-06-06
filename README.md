# lncRNA Analysis Pipeline
LncRNA framework for identification of lncRNA dysregulations from RNAseq data.

# Overview
The lncRNA framework pipeline provides a comprehensive solution for lncRNA analysis, encompassing preprocessing, alignment, transcript quantification, differential expression analysis, and functional annotation. This document outlines the steps for processing a single sample. For batch processing, refer to the scripts folder.
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
fastp -i SRR24709142_1.fastq.gz -I SRR24709142_2.fastq.gz -o SRR24709142_1_clean.fastq.gz -O SRR24709142_2_clean.fastq.gz -h SRR24709142_report.html --trim_poly_g --cut_tail --cut_mean_quality 20 --length_required 50
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

### Step 7a: Merge GTF
```bash
cat $(<mergelist.txt) > merged.gtf
```

### Step 8: Merge Transcripts
### Step 8a: Remove duplicate coordinate-based entries
```bash
awk '!seen[$1, $4, $5]++' merged.gtf > merged_without_duplicates.gtf
```
### Step 8b: Remove duplicate transcript ID-based entries
```bash
python Remove_Duplicate_transcript_ID.py merged_without_duplicates.gtf output_merged_without_duplicates.gtf
```
### Step 9: Remove Mitochondrial Chromosome
```bash
awk '/^>/{p=($0!="chrM MT")} p' human_genome/GRCh38.p14.genome.fa > human_genome/GRCh38.p14.genome_no_chrM.fa
```
### Step 10: Extract lncRNA sequences from genome fasta
```bash
gffread -w filtered_lncRNAs.fa -g human_genome/GRCh38.p14.genome_no_chrM.fa merged_without_duplicates.gtf
```

### Step 11: Assess Coding Potential
Download pre-built models and run:
```bash
cpat.py -g filtered_lncRNAs.fa -x Human_Hexamer.tsv -d Human_logitModel.RData -o cpat_output.txt
```
### Step 12: Cross-referencing with GTF file (filtered_gencode.gtf) to update annotations
The output file (cpat_output.txt) is filtered to retain only transcripts with low coding potential scores (i.e., non-coding). A typical threshold for non-coding classification is applied (e.g., coding potential score < 0.38 for human lncRNAs).
```bash
awk '$4 < 0.38' cpat_output.txt > high_confidence_lncRNAs.txt
```
This command filters the 4th column (coding potential score) and extracts only those classified as non-coding. The high-confidence lncRNAs are cross-referenced with GTF files (e.g., filtered_gencode.gtf) to update annotations, ensuring they align with transcript IDs of interest. The high-confidence lncRNAs identified in high_confidence_lncRNAs.txt are used as input for functional annotation (e.g., with FEELnc_filter.pl) to classify lncRNAs into various functional categories (e.g., intergenic or antisense).

### Step 13: Extract transcript IDs from the filtered CPAT output
```bash
awk '{print $1}' high_confidence_lncRNAs.txt > high_confidence_transcripts.txt
```
### Step 13_a: Cross-reference the filtered transcripts with the GTF file
```bash
grep -Ff high_confidence_transcripts.txt merged_without_duplicates.gtf > updated_filtered_gencode.gtf
```
### Step 13_b: Verify the updated GTF file (optional)
```bash
head updated_filtered_gencode.gtf
```

### Step 14: Quantify Transcript Abundance
```bash
featureCounts -a merged_without_duplicates.gtf -o SRR24709142_counts.txt -T 8 -p -B -C SRR24709142_Aligned.sortedByCoord.out.bam
```

### Step 15: Process Count Files
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

### Step 16: Differential Expression Analysis
Use R with DESeq2 for differential expression analysis and visualization. Refer to the R scripts in repository.

### Step 17: Heatmap
Refer to the R scripts in repository.

### Step 18: Functional Annotation
```bash
FEELnc_filter.pl -i updated_filtered_gencode.gtf -a gencode.v47.annotation.gtf -b transcript_biotype=protein_coding -l > candidate_lncRNA.gtf
awk 'BEGIN{FS="\t"; OFS="\t"} {match($9, /gene_id \"([^\"]+)\"/, a); if (a[1] != "") print a[1]}' candidate_lncRNA.gtf | sort | uniq > lncRNA_target_genes.txt
```

### Step 19: Pathway Enrichment Analysis
Use `clusterProfiler` in R for pathway enrichment analysis. Refer to the R scripts in repository.

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
