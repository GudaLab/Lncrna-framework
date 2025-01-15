#create a conda environment and install the below tools
FastQC (https://github.com/s-andrews/FastQC)
fastp (https://github.com/OpenGene/fastp) 
STAR (--quantMode GeneCounts) (https://github.com/alexdobin/STAR)
SAMtools (http://www.htslib.org/)
StringTie (#with -e -B for quantification)
StringTie Merge (-G for reference)
gffcompare (-G) (https://github.com/gpertea/gffcompare)
awk (keep u and x, length > 200 bp, exons ≥ 2)
CPAT (https://pypi.org/project/CPAT/)
Quantification: featureCounts (-p -B -C) (https://subread.sourceforge.net/featureCounts.html)
Differential Expression: DESeq2 (https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
Functional Annotation: FEELnc (https://github.com/tderrien/FEELnc)
Target Prediction: Miranda (https://github.com/hacktrackgnulinux/miranda)
Pathway Analysis: ClusterProfiler (enrichKEGG) (https://github.com/YuLab-SMU/clusterProfiler)
Visualization: ggplot2, igraph (for networks)
GENCODE Annotations for GTF/GFF Files (Download from here https://www.gencodegenes.org/human/)
Use GRCh38.p14 reference genome 

#Step 1: Download from SRA
fastq-dump SRR24709142 --split-3 --gzip

#Step 2: Quality Control of Raw Reads
fastqc SRR24709142_1.fastq.gz SRR24709142_2.fastq.gz -o qc_results

#Step 3: Trimming and Filtering (using fastp)
fastp -i SRR24709142_1.fastq.gz -I SRR24709142_2.fastq.gz -o SRR24709142_1_clean.fastq.gz -O SRR24709142_2_clean.fastq.gz -h SRR24709142_report.html --trim_poly_g --cut_mean_quality 20 --length_required 50

#Step 4: Star index
STAR --runThreadN 50 --runMode genomeGenerate --genomeDir human_genome/STAR_index --genomeFastaFiles human_genome/GRCh38.p14.genome.fa --sjdbGTFfile human_genome/gencode.v47.long_noncoding_RNAs.gtf

#Step 5: Alignment to Reference Genome
STAR --runThreadN 50 --genomeDir human_genome/STAR_index --readFilesIn SRR24709142_1_clean.fastq.gz SRR24709142_2_clean.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix SRR24709142_ --quantMode TranscriptomeSAM GeneCounts --outSAMstrandField intronMotif

#Step 6: Assemble Transcripts
cufflinks -g human_genome/gencode.v47.long_noncoding_RNAs.gtf -o SRR24709142_cufflinks_out -p 8 SRR24709142_Aligned.sortedByCoord.out.bam

#Step 7: Manually create a mergelist.txt file
vi mergelist.txt
press insert button to edit the file. 

#Type the following lines
SRR24709142_cufflinks_out/transcripts.gtf
SRR24709143_cufflinks_out/transcripts.gtf
SRR24709144_cufflinks_out/transcripts.gtf
SRR24709145_cufflinks_out/transcripts.gtf
SRR24709146_cufflinks_out/transcripts.gtf
SRR24709147_cufflinks_out/transcripts.gtf
SRR24709148_cufflinks_out/transcripts.gtf
SRR24709149_cufflinks_out/transcripts.gtf
SRR24709150_cufflinks_out/transcripts.gtf
SRR24709151_cufflinks_out/transcripts.gtf
SRR24709152_cufflinks_out/transcripts.gtf
SRR24709153_cufflinks_out/transcripts.gtf
SRR24709154_cufflinks_out/transcripts.gtf
SRR24709155_cufflinks_out/transcripts.gtf
SRR24709156_cufflinks_out/transcripts.gtf
SRR24709157_cufflinks_out/transcripts.gtf
SRR24709158_cufflinks_out/transcripts.gtf
SRR24709159_cufflinks_out/transcripts.gtf
SRR24709160_cufflinks_out/transcripts.gtf
SRR24709161_cufflinks_out/transcripts.gtf
SRR24709162_cufflinks_out/transcripts.gtf
SRR24709163_cufflinks_out/transcripts.gtf
SRR24709164_cufflinks_out/transcripts.gtf
SRR24709165_cufflinks_out/transcripts.gtf
SRR24709166_cufflinks_out/transcripts.gtf
SRR24709167_cufflinks_out/transcripts.gtf
SRR24709172_cufflinks_out/transcripts.gtf
SRR24709173_cufflinks_out/transcripts.gtf
SRR24709174_cufflinks_out/transcripts.gtf
SRR24709175_cufflinks_out/transcripts.gtf
SRR24709182_cufflinks_out/transcripts.gtf
SRR24709183_cufflinks_out/transcripts.gtf
SRR24709194_cufflinks_out/transcripts.gtf
SRR24709195_cufflinks_out/transcripts.gtf
SRR24709196_cufflinks_out/transcripts.gtf
SRR24709197_cufflinks_out/transcripts.gtf
SRR24709198_cufflinks_out/transcripts.gtf
SRR24709199_cufflinks_out/transcripts.gtf
SRR24709204_cufflinks_out/transcripts.gtf
SRR24709205_cufflinks_out/transcripts.gtf
SRR24709206_cufflinks_out/transcripts.gtf
SRR24709207_cufflinks_out/transcripts.gtf
SRR24709208_cufflinks_out/transcripts.gtf
SRR24709209_cufflinks_out/transcripts.gtf
SRR24709210_cufflinks_out/transcripts.gtf
SRR24709211_cufflinks_out/transcripts.gtf
SRR24709212_cufflinks_out/transcripts.gtf
SRR24709213_cufflinks_out/transcripts.gtf
SRR24709214_cufflinks_out/transcripts.gtf
SRR24709215_cufflinks_out/transcripts.gtf
SRR24709218_cufflinks_out/transcripts.gtf
SRR24709219_cufflinks_out/transcripts.gtf
SRR24709220_cufflinks_out/transcripts.gtf
SRR24709221_cufflinks_out/transcripts.gtf
SRR24709222_cufflinks_out/transcripts.gtf
SRR24709223_cufflinks_out/transcripts.gtf
SRR24709224_cufflinks_out/transcripts.gtf
SRR24709225_cufflinks_out/transcripts.gtf

press escape button and type :wq to save the file 

#Step 8: After gtf merge, clear duplicates
awk '!seen[$1, $4, $5]++' merged.gtf > merged_without_duplicates.gtf

#Step 9: Remove chrM from genome
awk '/^>/{p=($0!=">chrM MT")} p' human_genome/GRCh38.p14.genome.fa > human_genome/GRCh38.p14.genome_no_chrM.fa

#Step 10: Extract lncRNA sequences from genome fasta
gffread -w filtered_lncRNAs.fa -g human_genome/GRCh38.p14.genome_no_chrM.fa merged_without_duplicates.gtf

#Step 11: Assess Coding Potential (Download the pre-built models from here https://github.com/liguowang/cpat/tree/master/prebuilt_models)
cpat.py -g filtered_lncRNAs.fa -x Human_Hexamer.tsv -d Human_logitModel.RData -o cpat_output.txt

#Step 12: Cross-referencing with GTF file (filtered_gencode.gtf) to update annotations
awk '$4 < 0.38' cpat_output.txt > high_confidence_lncRNAs.txt

#Step 13: Extract transcript IDs from the filtered CPAT output
awk '{print $1}' high_confidence_lncRNAs.txt > high_confidence_transcripts.txt
#Step 13_a: Cross-reference the filtered transcripts with the GTF file
grep -Ff high_confidence_transcripts.txt merged_without_duplicates.gtf > updated_filtered_gencode.gtf
#Step 13_b: Verify the updated GTF file (optional)
head updated_filtered_gencode.gtf

#Step 14: Quantify Transcript Abundance
featureCounts -a merged_without_duplicates.gtf -o SRR24709142_counts.txt -T 8 -p -B -C SRR24709142_Aligned.sortedByCoord.out.bam

#Step 15:extract 1st and 7th column
for file in *_counts.txt; do
sed '1d' "$file" | awk -v OFS="\t" '{print $1, $7}' | sed 's/_Aligned\.sortedByCoord\.out\.bam//g' > "processed_$file"
done
#Extract the 1st row
cut -f1 processed_SRR24709142_counts.txt > combined_file.txt
#Loop through all other files and extract only the 2nd column
for file in processed_*_counts.txt; do
    awk '{print $2}' "$file" | paste combined_file.txt - > temp.txt
    mv temp.txt combined_file.txt
done

#Step 16: Differential Expression Analysis
setwd("C:/Users/s.jagadesan/Downloads")
library(DESeq2)
library(ggplot2)
library(ggrepel)
samples <- read.table("Luminal_A_vs_control.txt", header = TRUE)
count_data <- read.table("combined_file.txt", header = TRUE, row.names = 1)
samples <- samples[match(colnames(count_data), samples$sample_name), ]
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = samples, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
res$gene <- rownames(res)
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange >= 1 & res$padj <= 0.05] <- "UP"
res$diffexpressed[res$log2FoldChange <= -1 & res$padj <= 0.05] <- "DOWN"
top_genes <- rbind(
  res[res$diffexpressed == "UP", ][order(res[res$diffexpressed == "UP", ]$padj), ][1:25, ],
  res[res$diffexpressed == "DOWN", ][order(res[res$diffexpressed == "DOWN", ]$padj), ][1:25, ]
)
write.csv(as.data.frame(res), file = "differential_expression_results_with_labels.csv")

ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = diffexpressed)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 30) +
  theme_minimal() +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "gray")) +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 p-value", color = "Regulation") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

#Step 17:Heatmap
library(pheatmap)
normalized_counts <- counts(dds, normalized = TRUE)
significant_genes <- res[!is.na(res$padj) & res$padj <= 0.05, ]
significant_gene_names <- rownames(significant_genes)
heatmap_data <- normalized_counts[significant_gene_names, ]
heatmap_data_log <- log2(heatmap_data + 0.0001)
samples$group <- factor(samples$group, levels = unique(samples$group))
ordered_columns <- order(samples$group)
heatmap_data_log <- heatmap_data_log[, ordered_columns]
annotation_col <- data.frame(Cohort = samples$group[ordered_columns])
rownames(annotation_col) <- colnames(heatmap_data_log)
unique_cohorts <- levels(samples$group)
annotation_colors <- list(Cohort = setNames(rainbow(length(unique_cohorts)),  unique_cohorts))
pheatmap(
  heatmap_data_log,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  show_rownames = FALSE,
  show_colnames = TRUE, 
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of All Significant Genes"
)

#Heatmap top 25 up and downregulated genes
top_genes_names <- top_genes$gene 
heatmap_top_genes_data <- normalized_counts[top_genes_names, ]
heatmap_top_genes_data_log <- log2(heatmap_top_genes_data + 0.0001)
heatmap_top_genes_data_log <- heatmap_top_genes_data_log[, ordered_columns]
annotation_col <- data.frame(Cohort = samples$group[ordered_columns])
rownames(annotation_col) <- colnames(heatmap_top_genes_data_log)
unique_cohorts <- levels(samples$group)
annotation_colors <- list(Cohort = setNames(rainbow(length(unique_cohorts)), unique_cohorts))
pheatmap(
  heatmap_top_genes_data_log,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,  # Cluster genes (rows)
  cluster_cols = FALSE, # Do not cluster samples (columns)
  show_rownames = TRUE, # Show gene names for top 50 genes
  show_colnames = TRUE, # Show sample names
  color = colorRampPalette(c("blue", "white", "red"))(50), # Color gradient
  main = "Heatmap of Top 25 Up and Downregulated Genes"
)

#Step 18: Functional Annotation
FEELnc_filter.pl -i updated_filtered_gencode.gtf -a gencode.v47.annotation.gtf -b transcript_biotype=protein_coding -l > candidate_lncRNA.gtf
#extract_gene_name
awk 'BEGIN{FS="\t"; OFS="\t"} {match($9, /gene_id "([^"]+)"/, a); if (a[1] != "") print a[1]}' candidate_lncRNA.gtf | sort | uniq > lncRNA_target_genes.txt

#Step 19:Pathway Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
target_genes <- read.table("lncRNA_target_genes.txt", header = FALSE, stringsAsFactors = FALSE)
gene_symbols <- as.character(target_genes$V1)
gene_vector_entrez <- bitr(gene_symbols, fromType = "ENSEMBL",  toType = "ENTREZID",  OrgDb = org.Hs.eg.db)
enrich_res <- enrichKEGG(gene = gene_vector_entrez$ENTREZID,  organism = "hsa")
write.csv(as.data.frame(enrich_res@result), file = "pathway_enrichment_results_ENSEMBL_ID.csv")