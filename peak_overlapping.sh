#Cut&Run peaks overlap with RNAseq expression data

awk -F, 'NR==FNR{patterns[$0]=1; next} FNR==1{print; next} { temp=$NF; gsub(/"/, "", temp); if (patterns[temp]) print }' /fs/ess/PAS1348/nagesh/cutandrun_hannahs/output/mapping/bam/rmdups/macs3/H82_cMYC_25_75_ids.txt H82/Significant_DE_genes_H82_vs_H1048_batch_corrected.csv > H82/H82_gene_expression_associated_with_cMYC_peaks.csv
