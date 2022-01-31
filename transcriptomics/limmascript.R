setwd("/home/cristal/Documents/GEO")
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# To run this case study, you should have R version of 3.0.2 or later bioinf.wehi.edu.au/RNAseqCaseStudy/

# load libraries
library(Rsubread)
library(limma)
library(edgeR)

# read in target file
options(digits=2)
targets <- readTargets()
targets

# create a design matrix
celltype <- factor(targets$CellType)
design <- model.matrix(~celltype)

# build an index for reference sequence (Chr1 in hg19)
buildindex(basename="chr1",reference="hg19_chr1.fa")

# align reads
align(index="chr1",readfile1=targets$InputFile,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,unique=TRUE,indels=5)

# count numbers of reads mapped to NCBI Refseq genes
fc <- featureCounts(files=targets$OutputFile,annot.inbuilt="hg19")
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# generate RPKM values if you need them
x_rpkm <- rpkm(x,x$genes$Length)

# filter out low-count genes
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

# perform voom normalization
y <- voom(x,design,plot=TRUE)

# cluster libraries
plotMDS(y,xlim=c(-2.5,2.5))

# fit linear model and assess differential expression
fit <- eBayes(lmFit(y,design))
topTable(fit,coef=2)

#no Limma users guide, capítulo 15.4 diz
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

#Or, to give more weight to fold-changes in the gene ranking, one might use:
fit <- lmFit(logCPM, design)
fit <- treat(fit, lfc=log2(1.2))
topTreat(fit, coef=ncol(design))

# For Differential splicing, see chapter 15.7



  