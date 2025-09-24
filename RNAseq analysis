library(limma)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(EnhancedVolcano)

# Set publication-ready theme
theme_publication <- function() {
    theme_bw() +
    theme(
        # Text elements
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        legend.margin = margin(t = 15),
        
        # Panel and grid
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", color = NA),
        
        # Strip for facets
        strip.text = element_text(size = 12, face = "bold", margin = margin(5, 5, 5, 5)),
        strip.background = element_rect(fill = "grey95", color = "black", size = 1),
        
        # Margins
        plot.margin = margin(20, 20, 20, 20)
    )
}

# Extract gene IDs and expression matrix
gene_ids <- TPM_genes$gene_id
expr_data <- as.matrix(TPM_genes[, -1])
rownames(expr_data) <- gene_ids

# Clean gene symbol mapping - remove TXNAME and take unique rows
# This avoids duplicate gene symbols from multiple transcript variants
gene_mapping <- geneid_genesymbol %>%
    select(GENEID, GENESYMBOL) %>%
    distinct() %>%
    # Remove any NA or empty gene symbols
    filter(!is.na(GENESYMBOL), GENESYMBOL != "", GENESYMBOL != "NA") %>%
    # Sort by gene symbol for consistent processing
    arrange(GENESYMBOL, GENEID)

# Check for any remaining duplicates and keep first occurrence
# This handles cases where one gene ID maps to multiple symbols
gene_mapping <- gene_mapping %>%
    group_by(GENEID) %>%
    slice(1) %>%
    ungroup()

# Handle duplicate gene symbols (multiple gene IDs with same symbol)
# Add serial numbers to duplicate gene symbols
gene_mapping <- gene_mapping %>%
    group_by(GENESYMBOL) %>%
    mutate(
        gene_symbol_unique = if(n() > 1) {
            paste0(GENESYMBOL, "_", row_number())
        } else {
            GENESYMBOL
        }
    ) %>%
    ungroup()

# Create a clean lookup for gene symbols
gene_id_to_symbol <- setNames(gene_mapping$gene_symbol_unique, gene_mapping$GENEID)

# Map gene IDs in expression data to gene symbols
gene_symbols <- gene_id_to_symbol[rownames(expr_data)]

# For genes without mapping, keep the original gene ID
unmapped_genes <- is.na(gene_symbols)
gene_symbols[unmapped_genes] <- rownames(expr_data)[unmapped_genes]

# Ensure all gene symbols are unique using make.unique as final safety check
# This handles any edge cases where duplicates might still exist
gene_symbols <- make.unique(gene_symbols)

# Set the row names to gene symbols
rownames(expr_data) <- gene_symbols

# Create a clean mapping table for later use
# Ensure the data frame is properly constructed with unique row names
gene_symbol_to_id <- data.frame(
    gene_symbol = gene_symbols,
    gene_id = gene_ids,
    stringsAsFactors = FALSE
)

# Since gene_symbols are now guaranteed unique, this should work
rownames(gene_symbol_to_id) <- gene_symbols

cat("=== GENE SYMBOL MAPPING SUMMARY ===\n")
cat("Total gene IDs in expression data:", length(gene_ids), "\n")
cat("Total mappings in geneid_genesymbol (after cleaning):", nrow(gene_mapping), "\n")
cat("Successfully mapped to gene symbols:", sum(!unmapped_genes), "\n")
cat("Unmapped (kept as gene IDs):", sum(unmapped_genes), "\n")
cat("Duplicate gene symbols handled:", sum(gene_mapping$GENESYMBOL != gene_mapping$gene_symbol_unique), "\n")
cat("Final unique gene symbols:", length(unique(gene_symbols)), "\n\n")

# Define experimental groups
group <- factor(c(
    "H1048", "H1048", "H1048", "H1048", "H1048", "H1048", "H69", "H69", "H69", "H69"))

# Define batch information
batch <- factor(c(1, 1, 1, 2, 2, 2, 1, 2, 1, 1))

# Define sample names
sample_names <- c("H1048_1", "H1048_2", "H1048_3", "H1048_4", "H1048_5", "H1048_6",
                  "H69_1", "H69_2", "H69_3", "H69_4")
colnames(expr_data) <- sample_names

# Create design matrix
design <- model.matrix(~group)
colnames(design) <- c("Intercept", "groupH69")

# Filter genes with low expression by group
h1048_samples <- which(group == "H1048")
h69_samples <- which(group == "H69")
keep_h1048 <- rowSums(expr_data[, h1048_samples] > 1) >= 3
keep_h69 <- rowSums(expr_data[, h69_samples] > 1) >= 3
keep <- keep_h1048 | keep_h69
filtered_expr <- expr_data[keep, ]

cat("=== GENE FILTERING SUMMARY ===\n")
cat("Total genes before filtering:", nrow(expr_data), "\n")
cat("Genes kept after H1048 filter:", sum(keep_h1048), "\n")
cat("Genes kept after H69 filter:", sum(keep_h69), "\n")
cat("Genes kept after combined filter:", sum(keep), "\n\n")

# Create DGEList object and normalize
dge <- DGEList(counts = filtered_expr)
dge <- calcNormFactors(dge)

# Apply voom transformation
v <- voom(dge, design)

# Create batch-corrected data
cat("=== APPLYING BATCH EFFECT CORRECTION ===\n")
v_corrected <- v
v_corrected$E <- removeBatchEffect(v$E, batch = batch, design = design)
cat("Batch effect correction applied successfully\n\n")

# DIFFERENTIAL EXPRESSION ANALYSIS ON BATCH-CORRECTED DATA
cat("=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
cat("Performing DE analysis on batch-corrected data\n")

# Fit linear model on batch-corrected data
fit <- lmFit(v_corrected, design)

# Empirical Bayes moderation
fit <- eBayes(fit)

# Get results for H69 vs H1048 comparison
results_H69_vs_H1048 <- topTable(fit, coef = "groupH69", n = Inf, sort.by = "P")

# Add both gene symbols and gene IDs to results
results_H69_vs_H1048$gene_symbol <- rownames(results_H69_vs_H1048)
# Map back to original gene IDs for reference
symbol_to_id <- setNames(gene_mapping$GENEID, gene_mapping$GENESYMBOL)
results_H69_vs_H1048$gene_id <- symbol_to_id[results_H69_vs_H1048$gene_symbol]
# If no mapping found, use the symbol as ID (for unmapped genes)
results_H69_vs_H1048$gene_id[is.na(results_H69_vs_H1048$gene_id)] <- results_H69_vs_H1048$gene_symbol[is.na(results_H69_vs_H1048$gene_id)]

# Save full results
write.csv(results_H69_vs_H1048, "DE_results_H69_vs_H1048_batch_corrected.csv", row.names = FALSE)

# Extract significant genes (adjusted p-value < 0.05 and |logFC| > 2)
sig_genes <- results_H69_vs_H1048 %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 2) %>%
    arrange(P.Value)

# Ensure sig_genes has all the same columns as results_H69_vs_H1048
# including display_label
sig_genes <- sig_genes %>%
    mutate(
        display_label = ifelse(grepl("^ENSG", gene_symbol), gene_id, gene_symbol)
    )

# Save significant genes
write.csv(sig_genes, "Significant_DE_genes_H69_vs_H1048_batch_corrected.csv", row.names = FALSE)

# Note: Focusing on genes with |logFC| > 2 and adj.P.Val < 0.05

# VISUALIZATION OF RESULTS

# 1. Create Enhanced Volcano Plot
cat("Creating Enhanced Volcano Plot...\n")

# Prepare data for EnhancedVolcano
volcano_results <- results_H69_vs_H1048
rownames(volcano_results) <- volcano_results$gene_symbol

# Create enhanced volcano plot with all gene symbols
cat("Creating Enhanced Volcano Plot with gene symbols...\n")

# Prepare data for EnhancedVolcano - ensure row names match
volcano_results <- results_H69_vs_H1048
rownames(volcano_results) <- volcano_results$gene_symbol

# Filter significant genes and make sure they exist in volcano_results
sig_genes_labels <- sig_genes$gene_symbol
sig_genes_labels <- sig_genes_labels[sig_genes_labels %in% rownames(volcano_results)]

# Use display_label for better gene symbol representation
enhanced_volcano <- EnhancedVolcano(volcano_results,
    lab = rownames(volcano_results),  # Use rownames directly
    x = 'logFC',
    y = 'P.Value',
    title = 'H69 vs H1048 (Batch Corrected)',
    subtitle = 'Differential Expression Analysis',
    caption = paste0('Total genes = ', nrow(volcano_results)),
    pCutoff = 0.05,
    FCcutoff = 2.0,
    pointSize = 2.0,
    labSize = 3.0,
    labCol = 'black',
    colAlpha = 0.6,
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    max.overlaps = 50,  # Increased to show more gene symbols
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    selectLab = sig_genes_labels) +  # Label all significant genes
    theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10)
    )

# Save enhanced volcano plot
ggsave("Enhanced_Volcano_plot_H69_vs_H1048_batch_corrected.png", enhanced_volcano, 
       width = 12, height = 10, dpi = 300)

# 2. Create customized Enhanced Volcano with specific gene highlighting
# Get top 15 up and down regulated genes for highlighting
top_up <- sig_genes %>% filter(logFC > 0) %>% head(15)
top_down <- sig_genes %>% filter(logFC < 0) %>% head(15)

# Make sure the genes exist in volcano_results
highlight_genes <- c(top_up$gene_symbol, top_down$gene_symbol)
highlight_genes <- highlight_genes[highlight_genes %in% rownames(volcano_results)]

# Create custom volcano with selective labeling using proper gene symbols
custom_volcano <- EnhancedVolcano(volcano_results,
    lab = rownames(volcano_results),  # Use rownames directly for consistency
    x = 'logFC',
    y = 'P.Value',
    title = 'H69 vs H1048 (Batch Corrected) - Top DE Genes',
    subtitle = paste('Up-regulated:', nrow(top_up), '| Down-regulated:', nrow(top_down)),
    caption = paste0('FDR < 0.05, |log2FC| > 2'),
    pCutoff = 0.05,
    FCcutoff = 2.0,
    pointSize = 2.5,
    labSize = 3.5,
    selectLab = highlight_genes,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    colAlpha = 0.7,
    legendLabels = c('Not Sig.', 'Sig. FC', 'Sig. p-value', 'Sig. p-value & FC'),
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    colConnectors = 'black',
    max.overlaps = 30) +
    theme_publication()

ggsave("Custom_Enhanced_Volcano_H69_vs_H1048_batch_corrected.png", custom_volcano, 
       width = 14, height = 10, dpi = 300)

# 3. Create MA plot with updated fold change lines
ma_plot <- ggplot(results_H69_vs_H1048, aes(x = AveExpr, y = logFC)) +
    geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"),
                       labels = c("Not Significant", "Significant")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "blue", alpha = 0.7) +
    theme_publication() +
    labs(x = "Average log Expression", 
         y = "log2 Fold Change",
         title = "MA Plot: H69 vs H1048 (Batch Corrected)",
         color = "Significance") +
    theme(legend.position = "right")

ggsave("MA_plot_H69_vs_H1048_batch_corrected.png", ma_plot, 
       width = 10, height = 8, dpi = 300)

# 4. Create PCA plots before and after batch correction
pca_before <- prcomp(t(v$E))
pca_after <- prcomp(t(v_corrected$E))

# Calculate variance explained
var_before <- summary(pca_before)$importance[2, 1:2] * 100
var_after <- summary(pca_after)$importance[2, 1:2] * 100

# Create data frames for PCA plotting
pca_before_df <- data.frame(
    PC1 = pca_before$x[,1], 
    PC2 = pca_before$x[,2],
    Group = group,
    Batch = paste("Batch", batch),
    Sample = sample_names
)

pca_after_df <- data.frame(
    PC1 = pca_after$x[,1], 
    PC2 = pca_after$x[,2],
    Group = group,
    Batch = paste("Batch", batch),
    Sample = sample_names
)

# Define colors and shapes
group_colors <- c("H1048" = "#E31A1C", "H69" = "#1F78B4")
batch_shapes <- c("Batch 1" = 16, "Batch 2" = 17)

# Combined PCA plot (before and after)
pca_combined_df <- rbind(
    cbind(pca_before_df, Correction = "Before"),
    cbind(pca_after_df, Correction = "After")
)

combined_pca_plot <- ggplot(pca_combined_df, aes(x = PC1, y = PC2, color = Group, shape = Batch)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(aes(label = Sample), 
              nudge_x = 0, nudge_y = 1.5, 
              size = 2.5, color = "black", show.legend = FALSE) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = batch_shapes) +
    facet_wrap(~Correction, scales = "free") +
    theme_publication() +
    labs(
        x = paste0("PC1 (", round(var_before[1], 1), "% / ", round(var_after[1], 1), "%)"),
        y = paste0("PC2 (", round(var_before[2], 1), "% / ", round(var_after[2], 1), "%)"),
        title = "PCA Analysis: Before and After Batch Correction",
        color = "Cell Line",
        shape = "Batch"
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 4)),
        shape = guide_legend(override.aes = list(size = 4))
    )

ggsave("PCA_before_after_batch_correction_H69_vs_H1048.png", combined_pca_plot, 
       width = 14, height = 7, dpi = 300)

# SUMMARY STATISTICS
cat("=== DIFFERENTIAL EXPRESSION RESULTS SUMMARY ===\n")
cat("Total genes analyzed:", nrow(results_H69_vs_H1048), "\n")
cat("Significantly DE genes (adj.P.Val < 0.05, |logFC| > 2):", nrow(sig_genes), "\n")
cat("Up-regulated genes (logFC > 2):", sum(sig_genes$logFC > 2), "\n")
cat("Down-regulated genes (logFC < -2):", sum(sig_genes$logFC < -2), "\n")

# Top 10 up-regulated and down-regulated genes
cat("\nTop 10 Up-regulated genes (logFC > 2):\n")
top_up <- sig_genes %>% filter(logFC > 2) %>% head(10)
if(nrow(top_up) > 0) {
    print(top_up[c("gene_symbol", "gene_id", "logFC", "adj.P.Val")])
}

cat("\nTop 10 Down-regulated genes (logFC < -2):\n")
top_down <- sig_genes %>% filter(logFC < -2) %>% head(10)
if(nrow(top_down) > 0) {
    print(top_down[c("gene_symbol", "gene_id", "logFC", "adj.P.Val")])
}

cat("\n=== BATCH INFORMATION ===\n")
cat("Batch assignment:\n")
for(i in 1:length(sample_names)) {
    cat(sprintf("%s: %s, Batch %s\n", sample_names[i], group[i], batch[i]))
}

cat("\nBatch distribution:\n")
batch_dist <- table(group, batch)
print(batch_dist)

cat("\n=== FILES GENERATED ===\n")
cat("1. DE_results_H69_vs_H1048_batch_corrected.csv - Full DE results\n")
cat("2. Significant_DE_genes_H69_vs_H1048_batch_corrected.csv - Significant genes (|FC| > 2, p.adj < 0.05)\n")
cat("3. Enhanced_Volcano_plot_H69_vs_H1048_batch_corrected.png - Enhanced volcano plot with gene symbols\n")
cat("4. Custom_Enhanced_Volcano_H69_vs_H1048_batch_corrected.png - Custom volcano with top DE genes labeled\n")
cat("5. MA_plot_H69_vs_H1048_batch_corrected.png - MA plot\n")
cat("6. PCA_before_after_batch_correction_H69_vs_H1048.png - Combined PCA plot\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Differential expression analysis completed on batch-corrected data.\n")
cat("H69 vs H1048 comparison with |log2FC| > 2 threshold.\n")
cat("All results and visualizations have been saved.\n")
