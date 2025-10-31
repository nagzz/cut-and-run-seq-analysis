#!/usr/bin/env Rscript

# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(dplyr)
library(RColorBrewer)

# Create output directory with full path
output_dir <- file.path(getwd(), "combined_analysis_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check if directory was created successfully
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}
cat("Output directory created:", output_dir, "\n")

# Set working directory and file paths
setwd("/fs/ess/PAS1348/nagesh/novogene040725/output/H82")

# Define file paths
rnaseq_file <- "H82_all_gene_expression_associated_with_cMYC_peaks.csv"
bed_file <- "/fs/ess/PAS1348/nagesh/cutandrun_hannahs/output/mapping/bam/rmdups/macs3/H446_cMYC_transcripts_with_peaks_500_range_1.bed"
bigwig_file <- "/fs/ess/PAS1348/nagesh/cutandrun_hannahs/output/mapping/bam/rmdups/macs3/H446_c_MYC_bigwigcompare2_2024_08_27_12_35_48.bw"

# Function to read BED file and extract gene names
read_bed_with_genes <- function(bed_file) {
  bed_data <- fread(bed_file, header = FALSE)
  colnames(bed_data) <- c("chr", "start", "end", "name", "score", "strand", "start2", "end2")
  
  # Extract gene names and remove "_1" suffix
  bed_data$gene_symbol <- gsub("_1$", "", bed_data$name)
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = bed_data$chr,
    ranges = IRanges(start = bed_data$start, end = bed_data$end),
    strand = bed_data$strand,
    gene_symbol = bed_data$gene_symbol
  )
  
  return(list(gr = gr, data = bed_data))
}

# Function to create signal matrix from bigwig
create_signal_matrix <- function(regions, bigwig_file, upstream = 2000, downstream = 2000, bin_size = 50) {
  # Load bigwig
  bw <- import(bigwig_file, format = "BigWig")
  
  # Extend regions
  extended_regions <- resize(regions, width = upstream + downstream + width(regions), fix = "center")
  
  # Create bins
  n_bins <- (upstream + downstream) / bin_size
  bins <- tile(extended_regions, n = n_bins)
  
  # Calculate mean signal in each bin
  signal_matrix <- matrix(NA, nrow = length(regions), ncol = n_bins)
  
  for (i in 1:length(regions)) {
    region_bins <- bins[[i]]
    overlaps <- findOverlaps(region_bins, bw)
    
    for (j in 1:length(region_bins)) {
      hits <- overlaps[queryHits(overlaps) == j]
      if (length(hits) > 0) {
        signal_values <- bw[subjectHits(hits)]$score
        signal_matrix[i, j] <- mean(signal_values, na.rm = TRUE)
      } else {
        signal_matrix[i, j] <- 0
      }
    }
  }
  
  # Replace NAs with 0
  signal_matrix[is.na(signal_matrix)] <- 0
  
  return(signal_matrix)
}

cat("Reading data files...\n")
bed_info <- read_bed_with_genes(bed_file)
bed_data <- bed_info$data
rnaseq_data <- fread(rnaseq_file)

# Match genes between BED and RNA-seq data
common_genes <- intersect(bed_data$gene_symbol, rnaseq_data$gene_symbol)

# Filter and order data
bed_filtered <- bed_data[bed_data$gene_symbol %in% common_genes, ]
rnaseq_filtered <- rnaseq_data[rnaseq_data$gene_symbol %in% common_genes, ]
rnaseq_ordered <- rnaseq_filtered[match(bed_filtered$gene_symbol, rnaseq_filtered$gene_symbol), ]

# Create GRanges for filtered regions
regions_filtered <- GRanges(
  seqnames = bed_filtered$chr,
  ranges = IRanges(start = bed_filtered$start, end = bed_filtered$end),
  strand = bed_filtered$strand,
  gene_symbol = bed_filtered$gene_symbol
)

cat("Creating signal matrix...\n")
signal_matrix <- create_signal_matrix(regions_filtered, bigwig_file)
rownames(signal_matrix) <- bed_filtered$gene_symbol

# Prepare RNA-seq data for heatmap and sort by log2FC
rnaseq_matrix <- as.matrix(rnaseq_ordered$logFC)
rownames(rnaseq_matrix) <- rnaseq_ordered$gene_symbol
colnames(rnaseq_matrix) <- "log2FC"

# Sort by log2FC (descending order)
sort_order <- order(rnaseq_matrix[,1], decreasing = TRUE)
rnaseq_matrix_sorted <- rnaseq_matrix[sort_order, , drop = FALSE]
signal_matrix_sorted <- signal_matrix[sort_order, ]

cat("Creating combined plot...\n")

# Create the combined plot
pdf_file <- file.path(output_dir, "combined_CutandRunseq_rnaseq_plot.pdf")
cat("Creating PDF:", pdf_file, "\n")
pdf(pdf_file, width = 12, height = 8)

# ChIP-seq heatmap with gene symbol annotation
gene_annotation <- rowAnnotation(
  genes = anno_text(rownames(signal_matrix_sorted), 
                   gp = gpar(fontsize = 6),
                   location = unit(1, "npc"),
                   just = "right"),
  width = unit(4, "cm")
)

h1 <- Heatmap(signal_matrix_sorted,
              name = "ChIP-seq",
              col = colorRamp2(seq(min(signal_matrix), max(signal_matrix), length = 9),
                              brewer.pal(9, "YlOrRd")),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              column_title = "cMYC ChIP-seq Signal")

# RNA-seq heatmap
max_fc <- max(abs(rnaseq_matrix_sorted), na.rm = TRUE)
h2 <- Heatmap(rnaseq_matrix_sorted,
              name = "log2FC",
              col = colorRamp2(c(-max_fc, 0, max_fc), c("blue", "white", "red")),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              column_title = "RNA-seq log2FC",
              width = unit(2, "cm"))

# Draw combined plot
draw(gene_annotation + h1 + h2, column_title = paste0("Combined Analysis (", nrow(signal_matrix), " genes)"))

dev.off()

# Save summary statistics
summary_stats <- data.frame(
  total_genes_bed = nrow(bed_data),
  total_genes_rnaseq = nrow(rnaseq_data),
  common_genes = length(common_genes),
  upregulated = sum(rnaseq_matrix > 1, na.rm = TRUE),
  downregulated = sum(rnaseq_matrix < -1, na.rm = TRUE),
  unchanged = sum(abs(rnaseq_matrix) <= 1, na.rm = TRUE)
)

write.csv(summary_stats, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)

# Save gene order and data (sorted by log2FC)
gene_order_data <- data.frame(
  order = 1:length(common_genes),
  gene_symbol = rownames(rnaseq_matrix_sorted),
  logFC = rnaseq_matrix_sorted[,1],
  adj_pval = rnaseq_ordered$adj.P.Val[sort_order]
)

write.csv(gene_order_data, file.path(output_dir, "gene_order_and_data.csv"), row.names = FALSE)

cat("Analysis completed!\n")
cat("Output files in:", output_dir, "\n")
print(summary_stats)
