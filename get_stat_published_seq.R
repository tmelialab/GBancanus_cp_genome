#install.packages("read.gb")
#Tisha Melia
#This script to gather statistics of cp genomes

setwd("~/research/GBancanus/script/")
rm(list=ls())

library(Biostrings)
library(read.gb)

#-----------------------------------------------------------------------
# 1. Read in the fasta file
#-----------------------------------------------------------------------
dataSeq <- readDNAStringSet("input/EU849490_Gbancanus.fasta")
alphabetFrequency(dataSeq)
# A     C     G     T M R W S Y K V H D B   N - + .
# [1,] 51092 31240 30354 52027 0 0 0 0 0 0 0 0 0 0 397 0 0 0
width(dataSeq)#165110


dataSeqOur <- readDNAStringSet("../assembly/novoplasty/Option_1_Novop_Gonban_ILM_to_CPGAVAS2Rfcd176533.fas")
alphabetFrequency(dataSeqOur)
# A     C     G     T M R W S Y K V H D B N - + .
# [1,] 55405 32801 31770 56557 0 0 0 0 0 0 0 0 0 0 0 0 0 0
width(dataSeqOur) #176533
# 176533-165110
# [1] 11423


#-----------------------------------------------------------------------
# 2. Read in the gb file
#-----------------------------------------------------------------------
#read in the gb records
dataGB <- read.gb("../gb/EU849490.1.gb")
temp <- dataGB$EU849490$FEATURES; length(temp)#180
temp <- temp[names(temp) == "gene"]; length(temp)#79
geneNamesGB <- unname(unlist(lapply(temp, function(x){
  x[2,2]
}))); length(geneNamesGB)#79

dataGBOur <- read.gb("../gb/PQ046881_g_bancangus.gb")
temp <- dataGBOur$Gbancanus01_chloro$FEATURES
temp <- temp[names(temp) == "gene"]; length(temp)#139
geneNamesGBOur <- unname(unlist(lapply(temp, function(x){
  x[2,2]
}))); length(geneNamesGBOur)#139, 



#read g.affinis
dataGA <- read.gb("../gb/MN147872.1_G_affinis.gb")
temp <- dataGA$MN147872$FEATURES
temp <- temp[names(temp) == "gene"]; length(temp)#131
geneNamesGBOur <- unname(unlist(lapply(temp, function(x){
  x[2,2]
}))); length(geneNamesGBOur)#131, length 176548

176548-165110
131-79
79/131 = 0.603
