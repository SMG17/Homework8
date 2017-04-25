## Problem set 8

#2-a)
# Install DESeq2: (it estimates expression values, and calculates differential expression)
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
files <- c("C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output1/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output2/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output3/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output4/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output5/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output0/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output6/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output7/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output8/abundance.tsv",
           "C:/Users/SMG/AppData/Local/lxss/root/sequencing_class/HW8/output9/abundance.tsv")
names(files) <- c("sample1","sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10")

#2-b)
# Read in the Kallisto files using tximport:
txdat <- tximport(files, type = "kallisto", txOut = TRUE)

#2-c)
# Generate the condition matrix:
coldata <- data.frame(condition = c("SNF2","SNF2","SNF2","SNF2","SNF2","WT","WT","WT","WT","WT"))
rownames(coldata) = names(files)
coldata

#2-d)
# Turn it into a DESeq2 object:
dds <- DESeqDataSetFromTximport(txdat, colData = coldata, design =~ condition) #Tells it to group samples by the "condition" column of colData

# Run differential expression:
dds <- DESeq(dds)

# Summarize results:
res <- results(dds)

#3-a)
# A plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts. Points will be colored red if the adjusted p-value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
# An "MA-plot" provides a useful overview for an experiment with a two-group comparison (here the two conditions: WT and SNF2). The log2 fold change (using the mean of the counts for all replicates within WT, and within SNF2) for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor (by group/condition) is shown on the x-axis ("M" for minus, because a log ratio is equal to log minus log, and "A" for average).
# Each gene is represented with a dot. Genes with an adjusted p-value below a threshold (here 0.1, the default) are shown in red. The DESeq2 package incorporates a prior on log2 fold changes, resulting in moderated log2 fold changes from genes with low counts and highly variable counts, as can be seen by the narrowing of spread of points on the left side of the plot. This plot demonstrates that only genes with a large average normalized count contain sufficient information to yield a significant call.
plotMA(res) #, main="DESeq2", ylim=c(-2,2)

# Whether a gene is called significant depends not only on its LFC but also on its within-group variability, which DESeq2 quantifies as the dispersion. For strongly expressed genes, the dispersion can be understood as a squared coefficient of variation: a dispersion value of 0.01 means that the gene's expression tends to differ by typically ???0.01 = 10% between samples of the same treatment group. For weak genes, the Poisson noise is an additional source of noise.
# The function plotDispEsts visualizes DESeq2's dispersion estimates: each dot is a gene, and on the y-axis shows how much variation there is in how the gene is expressed in the group/condition.
# The black points are the dispersion estimates for each gene as obtained by considering the information from each gene separately. Unless one has many samples, these values fluctuate strongly around their true values. Therefore, we fit the red trend line, which shows the dispersions' dependence on the mean, and then shrink each gene's estimate towards the red line to obtain the final estimates (blue points) that are then used in the hypothesis test. The blue circles above the main "cloud" of points are genes which have high gene-wise dispersion estimates which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.
plotDispEsts(dds)

#3-b)
nom_pvalue_nb <- sum(res$pvalue < 0.05, na.rm = TRUE) # Should missing values (including NaN) be removed?
nom_pvalue_nb
# 2251

false_discovery_nb <- sum(res$padj < 0.05, na.rm = TRUE)
false_discovery_nb
# 1694