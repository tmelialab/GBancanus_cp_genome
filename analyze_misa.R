#Tisha Melia
#Get SSR statistics
setwd("~/research/GBancanus/misa/")
rm(list=ls())

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)   # for import() of GenBank
library(dplyr)
library(genbankr)



### 1. Load MISA SSR results ----
MISA_RESULT <- "gban.fasta.misa"   
stat <- read.table(MISA_RESULT, as.is = T, sep = "\t", header = T); dim(stat)#51 7
stat_range <- GRanges(
  seqnames = Rle(rep("Gonystylus bancanus", nrow(stat))),
  ranges = IRanges(stat$start, 
                   end = stat$end, 
                   names = stat$ID),
  strand = Rle(strand(rep("*", nrow(stat)))))


### 2. Define LSC / SSC / IR genomic partitions ----
annot_range <- GRanges(
  seqnames = Rle(rep("Gonystylus bancanus", 4)),
  ranges = IRanges(c(1, 88928, 131043, 134419), 
                   end = c(88927, 131042, 134418, 176533), 
                   names = c("LSC", "IRA", "SSC", "IRB")),
  strand = Rle(strand(rep("*", 4))))


countOverlaps(stat_range, annot_range, ignore.strand = T)
# LSC IRA SSC IRB 
# 35   8   0   8 
region_hits <- findOverlaps(stat_range, annot_range, ignore.strand = T)
region_hits <- data.frame(region_hits)
stat_range$region <- "LSC"
stat_range$region[region_hits$queryHits[region_hits$subjectHits == 2]] <- "IRA"
stat_range$region[region_hits$queryHits[region_hits$subjectHits == 3]] <- "SSC"
stat_range$region[region_hits$queryHits[region_hits$subjectHits == 4]] <- "IRB"
table(stat_range$region)
# IRA IRB LSC 
# 8   8  35 

### 3. Load GenBank annotation and extract genes ----
GB_FILE <- "../gb/PQ046881_g_bancangus.gb"
gb <- readGenBank(GB_FILE)
gb_coding <- gb@genes


### 4. Determine coding vs noncoding for each SSR ----
overlap_hits <- findOverlaps(stat_range, gb_coding, ignore.strand = T)
stat_range$coding_status <- "noncoding"
stat_range$coding_status[queryHits(overlap_hits)] <- "coding"

table(stat_range$coding_status)
# coding noncoding 
# 18        33 

### 5. Summary table ----
toprint <- data.frame(stat_range)
toprint$ID <- names(stat_range)
write.table(toprint, file = "ssr_properties.txt", sep = "\t")
