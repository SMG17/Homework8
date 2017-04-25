## Problem set 8

#2-a)
# Install DESeq2:
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
browseVignettes("DESeq2")

# Install GenomeInfoDb (independently of DESeq2):
source("https://bioconductor.org/biocLite.R")
biocLite("GenomeInfoDb")

# Install tximport:
source("https://bioconductor.org/biocLite.R")
biocLite("tximport")

# Install rhdf5 (independently of tximport):
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")

#2-b)
# Load libraries:
library(DESeq2)
library(GenomeInfoDb)
library(tximport)
library(rhdf5)

# Create a named vector of files to import:
files <- c("C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output0/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output1/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output2/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output3/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output4/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output5/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output6/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output7/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output8/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output9/abundance.tsv")
names(files) <- c("sample1","sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10")

# Read in the Kallisto files using tximport:
txdat <- tximport(files, type = "kallisto", txOut = TRUE)

# Generate the condition matrix:
coldata <- data.frame(condition = c("WT","SNF2","SNF2","SNF2","SNF2","SNF2","WT","WT","WT","WT"))
rownames(coldata) = names(files)
coldata

# Turn it into a DESeq2 object:
dds <- DESeqDataSetFromTximport(txdat, colData = coldata, design =~ condition) #Tells it to group samples by the "condition" column of colData

# Run differential expression:
dds <- DESeq(dds)

# Summarize results:
res <- results(dds)