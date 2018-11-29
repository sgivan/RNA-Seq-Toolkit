#!/bin/env Rscript
#SBATCH -J $prefix
#SBATCH -o ${prefix}_%J.o
#SBATCH -e ${prefix}_%J.e 
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=$memory

library(edgeR)
library(DESeq2)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(vegan)
library(pheatmap)

# Thanks Stephen Turner (see https://gist.github.com/stephenturner/4a599dbf120f380d38e7#file-volcanoplot-r)
plot_Volcano <- function(DE_results_filename,labelPoints=FALSE) {

    res <- read.table(DE_results_filename, na.strings="", sep="\t",header=TRUE)
    row.names(res) <- res$$Gene

    # No NA's in log2FoldChange, right?
    # True is also 1, right?
    log2fc <- res$$log2FoldChange
    circle_sizes <- rep(1,length(log2fc))

    res <- res%>%mutate(color = ifelse(abs(log2FoldChange) >= 1,ifelse(padj<=0.05 , "green", "orange"),ifelse(padj<=0.05,"red","black")))
    max_y <- max(-log10(res$$pvalue)) + 0.5

    ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue), ymax=max_y, ymin=0)) +
    geom_point(aes(colour=color), size=2.5) +
    scale_colour_manual(name="", values = c("green"="green", "orange"="orange", "red"="red", "black"="black")) +

    theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) + # Make axis text and title bigger
    guides(color=FALSE)

}
#
# The gene_cnt_matrix.tab file, loaded in the next command, should have gene IDs in the first column.
# Then, each sample should occupy its own column.
# Control coumns should be the columns after the first.
# Experimental columns should follow after all the control columns.
# the header should reflect the columns. For example:
# GeneID    Sample_1    Sample_2    Sample_3    Sample_4
# In this example, Sample_1 and Sample_2 are controls; Sample_3 and Sample_4 are the experimentals.
data = read.delim("$datafile", header=T, row.names=1, sep=" ", stringsAsFactors=F)
# may need round() for RSEM data. Not needed for STAR.
#rnaseqMatrix = round(data)
rnaseqMatrix = data

# Keep only rows with at least two cpms greater than one
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]

# Create a dataframe of condition labels
#condition = data.frame(condition=factor(c(rep("D", 5), rep("C",5))))
# $cntCont and $cntExp will be replaced by the number of Control samples and the number
# of Experimental samples, respectively.
condition = data.frame(condition=factor(c(rep("Cont", $cntCont), rep("Exp",$cntExp))))
sample_IDs <- colnames(rnaseqMatrix)

# Add rownames that correspond to sample IDs
rownames(condition) = sample_IDs

# DESeq2: Calculate differential expression using count data
ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = rnaseqMatrix,
    colData = condition,
    design = ~ condition)
dds = DESeq(ddsFullCountTable)

# plotPCA can't take dds straight
rldds <- rlog(dds)

pdf(paste0("$prefix", "_PCA.pdf"))
plotPCA(rldds,intgroup=c('condition')) # This will give groups, but not sample names within the group
dev.off()

svg(paste0("$prefix", "_PCA.svg"))
plotPCA(rldds,intgroup=c('condition')) # This will give groups, but not sample names within the group
dev.off()

#contrast=c("condition","D","C")
contrast=c("condition","Exp","Cont")
res = results(dds, contrast)
dds.res <- res

res = cbind(experimental="Exp", control="Cont", as.data.frame(res))

# Add baseMean Exp and Control
baseMeanExp     <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$$condition == "Exp"])
baseMeanControl <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$$condition == "Cont"])
res = cbind(baseMeanExp, baseMeanControl, as.data.frame(res))

res$$padj[is.na(res$$padj)]  <- 1
geneSymbol <- mapIds(org.Mm.eg.db, keys=row.names(res), column="SYMBOL",keytype="ENSEMBL", multiVals="first")
bestGeneDescriptor <- geneSymbol
bestGeneDescriptor[ is.na(bestGeneDescriptor) ] <- row.names(res)[ is.na(bestGeneDescriptor) ]

res <- cbind(GeneSymbol=geneSymbol, as.data.frame(res))
res <- cbind(Gene=row.names(res), as.data.frame(res))
res <- cbind(BestGeneDescriptor=bestGeneDescriptor, as.data.frame(res))

#write.table(as.data.frame(res[order(res$$pvalue),]), file='DESeq2_DE_results.txt', sep="\t", na="", row.names=FALSE,quote=FALSE)
write.table(as.data.frame(res[order(res$$pvalue),]), file=paste0("$prefix", "_DESeq2_DE_results.txt"), sep="\t", na="", row.names=FALSE,quote=FALSE)
write.table(rnaseqMatrix, file=paste0("$prefix", "_DESeq2_count_matrix.txt"), sep="\t", na="", row.names=FALSE,quote=FALSE)

pdf(paste0("$prefix", "_Volcano.pdf"))
plot_Volcano(paste0("$prefix", "_DESeq2_DE_results.txt"))
dev.off()

svg(paste0("$prefix", "_Volcano.svg"))
plot_Volcano(paste0("$prefix", "_DESeq2_DE_results.txt"))
dev.off()

vsd <- vst(dds, blind=F)
vsd.df.t <- t(as.data.frame(assay(vsd)))
vsd.adonis <- adonis(vsd.df.t ~ colData(dds)$$condition, method="eu", permutations=10000)

sink(paste0("$prefix","_adonis.txt"))
vsd.adonis
sink()

select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:1000]
pdf(paste0("$prefix", "_heatmap.pdf"))
pheatmap(assay(rldds)[select,], cluster_rows=T, cluster_cols=T, show_rownames=F, annotation_col=condition)
dev.off()

svg(paste0("$prefix", "_heatmap.svg"))
pheatmap(assay(rldds)[select,], cluster_rows=T, cluster_cols=T, show_rownames=F, annotation_col=condition)
dev.off()

save.image(file=paste0("$prefix", "_RData"))

