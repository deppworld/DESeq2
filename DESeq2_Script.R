#Sample condition N (normal), T (Treated)
#SRR2677524_N_htseq_count.tabular, SRR2677524_N_htseq_count[1].tabular, SRR2678264_T_htseq_count.tabular, SRR2678360_T_htseq_count.tabular, SRR2678360_T_htseq_count[1].tabular, SRR2678360_T_htseq_count[2].tabular, SRR2678360_T_htseq_count[3].tabular

#Set directory or you can skip this step if you are working in same directory, change the path according to your folder.
directory <- "C:/Users/dverma2/Desktop/HTSeq"
set(directory)
#Creating sample Matrix and input files from HTSeq data
sampleFiles <- grep("count",list.files(directory),value=TRUE)
sampleCondition = sapply(strsplit(sampleFiles, "_", fixed=T), function(x) x[2])#condition is separated by "." or "_"
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
#DE analysis
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition) #creating DESeq data set
ddsHTSeq
keep <- rowSums(counts(ddsHTSeq)) >= 10 #Optional step: a minimal pre-filtering to keep only rows that have at least 10 reads total
dds <- ddsHTSeq[keep,]

#command for DE
dds <- DESeq(dds)
res <- results(dds)
res                             #result file
#Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_T_vs_N", type="apeglm")
resLFC
#p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]
resOrdered
summary(res)  #To see the summary of pval and adj-pval
#MA-plot: plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2)) # MA-plot for the shrunken log2 fold changes, You can increase the FC limi like (10,-10), (15, -15)
# Data transformations and visualization
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
#Heatmap of the count matrix

library("pheatmap")
df <- as.data.frame(colData(dds))
pheatmap(assay(ntd), cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#Approach to count outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
#Dispersion plot and fitting alternatives
plotDispEsts(dds)
