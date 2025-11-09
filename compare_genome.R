#Tisha Melia
#12 October 2025
#This script compares the chloroplast genome: EU849490.1 and PQ046881
rm(list=ls())
setwd("/Users/tisha/research/GBancanus/script")

library(seqinr)
library(Biostrings)
library(genbankr)

# ---- INPUT FILES ----
old_fasta <- "../gb/EU849490.fasta"
new_fasta <- "../gb/PQ046881.fasta"
old_gb <- "../gb/EU849490.1.gb"
new_gb <- "../gb/PQ046881_g_bancangus.gb"

# ---- Function to summarize assembly stats ----
get_fasta_stats <- function(fasta_file) {
  seq <- readDNAStringSet(fasta_file)
  seq_length <- width(seq)
  total_length <- sum(seq_length)
  gc_count <- sum(letterFrequency(seq, c("G", "C")))
  gc_content <- round((gc_count / total_length) * 100, 2)
  n_count <- sum(letterFrequency(seq, "N"))
  
  return(data.frame(
    Length_bp = total_length,
    GC_content_percent = gc_content,
    N_count = n_count
  ))
}

# ---- Function to summarize annotation stats ----
get_gb_stats <- function(gb_file) {
  gb <- readGenBank(gb_file)
  # Extract gene-level info
  genes <- gb@genes
  
  # Extract gene names and IDs
  gene_names <- as.character(mcols(genes)$gene)
  gene_ids <- as.character(mcols(genes)$gene_id)
  
  # Normalize text for matching
  all_names <- tolower(paste(gene_names, gene_ids))
  
  # Count tRNAs and rRNAs by pattern
  trna_count <- sum(grepl("trn", all_names, ignore.case = TRUE))
  rrna_count <- sum(grepl("rrn|16s|23s|5s|4.5s", all_names, ignore.case = TRUE))

  # Everything else as protein-coding
  total_genes <- length(unique(all_names))
  protein_coding <- total_genes - (trna_count + rrna_count)
  
  return(data.frame(
    Annotated_genes = total_genes,
    Protein_coding_genes = protein_coding,
    tRNAs = trna_count,
    rRNAs = rrna_count
  ))
}

# ---- Run analyses ----
old_stats <- cbind(Source = "EU849490", get_fasta_stats(old_fasta), get_gb_stats(old_gb))
new_stats <- cbind(Source = "PQ046881", get_fasta_stats(new_fasta), get_gb_stats(new_gb))

comparison <- rbind(old_stats, new_stats)

# ---- Add improvement comments manually if desired ----
comparison$Comment <- c("Partial assembly", "Complete assembly with no gaps")

# ---- Save as CSV ----
write.csv(comparison, "output/chloroplast_assembly_comparison.csv", row.names = FALSE)

# ---- Print summary ----
print(comparison)


