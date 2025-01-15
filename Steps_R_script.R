#Step 14: Differential Expression Analysis
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
ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = diffexpressed)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, box.padding = 0.5, point.padding = 0.3, max.overlaps = 30) +
  theme_minimal() +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "gray")) +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 p-value", color = "Regulation") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
write.csv(as.data.frame(res), file = "differential_expression_results_with_labels.csv")

#Histogram of log2 Fold Changes
ggplot(as.data.frame(res), aes(x = log2FoldChange)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Histogram of Log2 Fold Changes",
    x = "Log2 Fold Change",
    y = "Frequency"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

#Step 15:Heatmap
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

#Step 16: Functional Annotation
FEELnc_filter.pl -i filtered_gencode.gtf -a gencode.v47.annotation.gtf -b transcript_biotype=protein_coding -l > candidate_lncRNA.gtf
#extract_gene_name
awk 'BEGIN{FS="\t"; OFS="\t"} {match($9, /gene_id "([^"]+)"/, a); if (a[1] != "") print a[1]}' candidate_lncRNA.gtf | sort | uniq > lncRNA_target_genes.txt

#Step 17:Pathway Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
target_genes <- read.table("lncRNA_target_genes.txt", header = FALSE, stringsAsFactors = FALSE)
gene_symbols <- as.character(target_genes$V1)
gene_vector_entrez <- bitr(gene_symbols, fromType = "ENSEMBL",  toType = "ENTREZID",  OrgDb = org.Hs.eg.db)
enrich_res <- enrichKEGG(gene = gene_vector_entrez$ENTREZID,  organism = "hsa")
write.csv(as.data.frame(enrich_res@result), file = "pathway_enrichment_results_ENSEMBL_ID.csv")