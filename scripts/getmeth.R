# Load necessary libraries
library(GenomicRanges)
library(genomation)
library(rtracklayer)
library(data.table)
library(GenomicFeatures)
library(dplyr)

# Define the paths to the GFF3 file and cov file
gff3_file <- "/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/Mus_musculus.GRCm38.96.gff3"
cov_file <- "/scratch/workspace/michael_olufemi_student_uml_edu-nanopore/space_omics/mapping/GLDS-48_M26_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz"

# Load the GFF3 file and create TxDb object
txdb <- makeTxDbFromGFF(gff3_file, format = "gff3")

# Get gene annotations
genes <- genes(txdb)

# Define promoter regions based on strand
# Promoters are 1500 bp upstream and 500 bp downstream for positive strand genes
# For negative strand genes, the upstream/downstream is reversed
promoters <- promoters(genes, upstream = 1500, downstream = 500)

# Adjust promoter regions based on strand orientation
promoters_positive <- promoters[strand(genes) == "+",]
promoters_positive <- promoters(promoters_positive, upstream = 1500, downstream = 500)  # For positive strand

promoters_negative <- promoters[strand(genes) == "-",]
promoters_negative <- promoters(promoters_negative, upstream = 500, downstream = 1500)  # For negative strand

# Combine promoter regions by strand
promoters <- c(promoters_positive, promoters_negative)

# Load the methylation cov file
methylation_summary <- fread(cov_file, col.names = c("chrom", "start", "end", "pct", "numCs", "numTs"))

# Function to calculate gene and promoter methylation summaries based on gene strand
calculate_gene_promoter_methylation <- function(methylation_summary, genes, promoters) {
  gr <- GRanges(seqnames = methylation_summary$chrom,
                ranges = IRanges(start = methylation_summary$start, end = methylation_summary$end),
                meth_5mC = methylation_summary$pct)  # Using the methylation percentage
  
  # Calculate gene body methylation
  gene_overlaps <- findOverlaps(gr, genes)
  gene_methylation <- data.frame(
    gene_id = names(genes)[subjectHits(gene_overlaps)],
    meth_5mC = gr$meth_5mC[queryHits(gene_overlaps)]
  )
  
  gene_methylation_summary <- gene_methylation %>%
    group_by(gene_id) %>%
    summarise(gene_body_methylation = mean(meth_5mC, na.rm = TRUE))
  
  # Calculate promoter methylation
  promoter_overlaps <- findOverlaps(gr, promoters)
  promoter_methylation <- data.frame(
    gene_id = names(promoters)[subjectHits(promoter_overlaps)],  # Match gene ID from promoters
    meth_5mC = gr$meth_5mC[queryHits(promoter_overlaps)]
  )
  
  promoter_methylation_summary <- promoter_methylation %>%
    group_by(gene_id) %>%
    summarise(promoter_methylation = mean(meth_5mC, na.rm = TRUE))
  
  # Merge gene body and promoter methylation summaries by gene ID
  combined_methylation <- merge(gene_methylation_summary, promoter_methylation_summary, by = "gene_id", all = TRUE)
  
  return(combined_methylation)
}

# Run the methylation calculation for the provided cov file
combined_methylation <- calculate_gene_promoter_methylation(methylation_summary, genes, promoters)

# Output the final table
write.csv(combined_methylation, "gene_promoter_body_methylation.csv", row.names = FALSE)

print("Gene promoter and gene body methylation summary created.")
