#!/usr/bin/env Rscript
#PBS -N $prefix
#PBS -d ./
#PBS -l nodes=1:ppn=1,mem=$memory
#PBS -q $jobqueue
#PBS -V

# cd PBS_O_WORKDIR

library(tidyverse)
library(edgeR)
library(DESeq2)
library($orgdb)
#library(ggplot2)
#library(dplyr)
library(vegan)
library(pheatmap)
library(gprofiler2)

# Thanks Stephen Turner (see https://gist.github.com/stephenturner/4a599dbf120f380d38e7#file-volcanoplot-r)
plot_Volcano <- function(DE_results_filename,labelPoints=FALSE) {

#    res <- read.table(DE_results_filename, na.strings="", sep="\t",header=TRUE)
    res <- read_tsv(DE_results_filename, col_names=T)
#    row.names(res) <- res$$Gene

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

# create the merged data object
merged.data <- read_tsv(paste0("../","$aligndir","/Sample_1/PE_ReadsPerGene.out.tab"), col_names=F, skip=4)[,c(1,$strand)]
sampnames <- c("GENEID", "Sample_1")
eod <- $cntCont + $cntExp
for (i in 2:eod) {
    sname <- paste0("Sample_",i)
    fpth <- paste0("../", "$aligndir","/", sname, "/PE_ReadsPerGene.out.tab")
#    sampdata <- read.delim(fpth, header=T, row.names=NULL, sep="\t", stringsAsFactors=F)
    sampnames <- append(sampnames, sname)
    sampdata <- read_tsv(fpth, col_names=F, skip=4)[,c(1,$strand)]
    merged.data <- merge(merged.data, sampdata, by="X1")
}
colnames(merged.data) <- sampnames
#
# The gene_cnt_matrix.tab file, loaded in the next command, should have gene IDs in the first column.
# Then, each sample should occupy its own column.
# Control coumns should be the columns after the first.
# Experimental columns should follow after all the control columns.
# the header should reflect the columns. For example:
# GeneID    Sample_1    Sample_2    Sample_3    Sample_4
# In this example, Sample_1 and Sample_2 are controls; Sample_3 and Sample_4 are the experimentals.
#data = read.delim("$datafile", header=T, row.names=1, sep=" ", stringsAsFactors=F)
data <- tibble::column_to_rownames(merged.data, var="GENEID")
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

#
# wrap the call to mapIds in a tryCatch() routine to catch errors
rslt <- tryCatch(
                    expr = {
                        geneSymbol <- mapIds($orgdb, keys=row.names(res), column="SYMBOL", keytype="$dbkey", multiVals="first")
                    },
                    error = function(e){
                        message("mapping gene symbols with org DB generated an error")
                        return(e)
                    },
                    warning = function(w){
                        message("mapping gene symbols with org DB generated a warning")
                        print(w)
                    },
                    finally = {
                        message("finished mapping gene symbols with org DB")
                    }
                 )

#
# if mapIds() generated an error, try to get gene names using the gprofiler2 function gconvert
if (inherits(rslt, 'simpleError')) {
    message("because of error using org DB, now try gprofiler2 to get gene symbols")
    gcrslt <- gconvert(row.names(res), organism="$gProfilerkey")
    geneSymbol <- gcrslt$$name
}

#
# let's assume that if we've got this far that geneSymbol contains gene names
bestGeneDescriptor <- gsub("'","_",geneSymbol)
bestGeneDescriptor[ is.na(bestGeneDescriptor) ] <- row.names(res)[ is.na(bestGeneDescriptor) ]

res <- cbind(GeneSymbol=geneSymbol, as.data.frame(res))
res <- cbind(Gene=row.names(res), as.data.frame(res))
res <- cbind(BestGeneDescriptor=bestGeneDescriptor, as.data.frame(res))

#
# write DE results and count matrix to separate files
write.table(as.data.frame(res[order(res$$pvalue),]), file=paste0("$prefix", "_DESeq2_DE_results.txt"), sep="\t", na="", row.names=FALSE,quote=FALSE)
write.table(rnaseqMatrix, file=paste0("$prefix", "_DESeq2_count_matrix.txt"), sep="\t", na="", row.names=FALSE,quote=FALSE)

#
# generate a volcano plot as a PDF file
pdf(paste0("$prefix", "_Volcano.pdf"))
plot_Volcano(paste0("$prefix", "_DESeq2_DE_results.txt"))
dev.off()

#
# generate a volcano plot as a SVG file
svg(paste0("$prefix", "_Volcano.svg"))
plot_Volcano(paste0("$prefix", "_DESeq2_DE_results.txt"))
dev.off()

#
# calcluate statistical significance of clustering in PCA analysis
vsd <- vst(dds, blind=F)
vsd.df.t <- t(as.data.frame(assay(vsd)))
vsd.adonis <- adonis(vsd.df.t ~ colData(dds)$$condition, method="eu", permutations=10000)

#
# output adonis() results to text file
sink(paste0("$prefix","_adonis.txt"))
vsd.adonis
sink()

# make heatmap plots

# pdf file
select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:1000]
pdf(paste0("$prefix", "_Heatmap.pdf"))
pheatmap(assay(rldds)[select,], cluster_rows=T, cluster_cols=T, show_rownames=F, annotation_col=condition)
dev.off()

# svg file
svg(paste0("$prefix", "_Heatmap.svg"))
pheatmap(assay(rldds)[select,], cluster_rows=T, cluster_cols=T, show_rownames=F, annotation_col=condition)
dev.off()

# do pathway analysis
#library(gProfileR)
#library(gprofiler2)
BGD.05 <- as.character(dplyr::select(dplyr::filter(dplyr::arrange(res, padj), padj < 0.05), BestGeneDescriptor)$$BestGeneDescriptor)
#gprofile_Ordered <- gprofiler(BGD.05, organism="$gProfilerkey", ordered_query=T, correction_method='analytical', sort_by_structure=T, significant=T)
gprofile_Ordered <- gost(BGD.05, organism="$gProfilerkey", ordered_query=T, correction_method='gSCS', significant=T)
#
# what we want to output contains a list, so collapse it
gprofile_Ordered$$result$$parents <- paste(gprofile_Ordered$$result$$parents, collapse=',')
# write the gprofiler output to a file
write.table(gprofile_Ordered['result'], file=paste0("$prefix","_gProfileR.txt"), sep="\t", quote=F, row.name=F, col.name=T)

#
# save image file
save.image(file=paste0("$prefix", "_RData"))

