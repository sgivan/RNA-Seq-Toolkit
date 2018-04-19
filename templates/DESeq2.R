#!/bin/env Rscript
#SBATCH -J DESeq2
#SBATCH -o out_DESeq2.o_%j
#SBATCH -e out_DESeq2.e_%j
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=40G

library(edgeR)
library(DESeq2)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)

# Thanks Stephen Turner (see https://gist.github.com/stephenturner/4a599dbf120f380d38e7#file-volcanoplot-r)
plot_Volcano <- function(DE_results_filename,labelPoints=FALSE) {

    res <- read.table(DE_results_filename, na.strings="", sep="\t",header=TRUE)
    row.names(res) <- res$Gene

    # No NA's in log2FoldChange, right?
    # True is also 1, right?
    log2fc <- res$log2FoldChange
    circle_sizes <- rep(1,length(log2fc))

    res <- res%>%mutate(color = ifelse(abs(log2FoldChange) >= 2,ifelse(padj<=0.05 , "green", "orange"),ifelse(padj<=0.05,"red","black")))
    max_y <- max(-log10(res$pvalue)) + 0.5

    ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue), ymax=max_y, ymin=0)) +
    geom_point(aes(colour=color), size=2.5) +
    scale_colour_manual(name="", values = c("green"="green", "orange"="orange", "red"="red", "black"="black")) +

    theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) + # Make axis text and title bigger
    guides(color=TRUE)

}

data = read.table("gene_cnt_matrix.tab", header=TRUE, row.names=1, com='', na.strings="", sep="\t")
#rnaseqMatrix = round(data)
rnaseqMatrix = data

# Keep only rows with at least two cpms greater than one
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]

# Create a dataframe of condition labels
#condition = data.frame(condition=factor(c(rep("D", 5), rep("C",5))))
condition = data.frame(condition=factor(c(rep("Cont", $cntCont), rep("Exp",$numExp))))
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

pdf("D_vs_C__DESeq2PCA_rldds.pdf")
plotPCA(rldds,intgroup=c('condition')) # This will give groups, but not sample names within the group
dev.off()

png("D_vs_C__DESeq2PCA_1200x1200.png",width=1200,height=1200)
plotPCA(rldds,intgroup=c('condition')) # This will give groups, but not sample names within the group
dev.off()

png("D_vs_C__DESeq2PCA_rldds.png",width=480,height=480)
plotPCA(rldds,intgroup=c('condition')) # This will give groups, but not sample names within the group
dev.off()

contrast=c("condition","D","C")
res = results(dds, contrast)

res = cbind(experimental="D", control="C", as.data.frame(res))

# Add baseMean Exp and Control
baseMeanExp     <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "D"])
baseMeanControl <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$condition == "C"])
res = cbind(baseMeanExp, baseMeanControl, as.data.frame(res))

res$padj[is.na(res$padj)]  <- 1
geneSymbol <- mapIds(org.Mm.eg.db, keys=row.names(res), column="SYMBOL",keytype="ENSEMBL", multiVals="first")
bestGeneDescriptor <- geneSymbol
bestGeneDescriptor[ is.na(bestGeneDescriptor) ] <- row.names(res)[ is.na(bestGeneDescriptor) ]

res <- cbind(GeneSymbol=geneSymbol, as.data.frame(res))
res <- cbind(Gene=row.names(res), as.data.frame(res))
res <- cbind(BestGeneDescriptor=bestGeneDescriptor, as.data.frame(res))

write.table(as.data.frame(res[order(res$pvalue),]), file='D_vs_C__DESeq2.DE_results', sep="\t", na="", row.names=FALSE,quote=FALSE)
write.table(rnaseqMatrix, file='D_vs_C__out.count_matrix', sep="\t", na="", row.names=FALSE,quote=FALSE)

pdf("D_vs_C__Volcano.pdf")
plot_Volcano('D_vs_C__DESeq2.DE_results')
dev.off()

png("D_vs_C__Volcano.png",width=480,height=480)
plot_Volcano('D_vs_C__DESeq2.DE_results')
dev.off()

png("D_vs_C__Volcano_1200x1200.png",width=1200,height=1200)
plot_Volcano('D_vs_C__DESeq2.DE_results')
dev.off()

