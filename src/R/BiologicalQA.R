#' ---
#' title: "Poplar (T89 - OP42) Cambium adventitious roots studies biological QA"
#' author: "Nicolas Delhomme & Alok Ranjan"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create palettes
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Process

#' ## All samples

#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/bellini_cambium/doc/samples.csv")

#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

# Rename the HTSeq files
names(res) <- samples$SampleName[match(sub("(_L[1,2])?_sortmerna.*\\.txt","",dir("htseq",pattern="*.txt")),samples$BGI_File_Name)]

#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
dir.create(file.path("analysis","HTSeq-new"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"analysis/HTSeq-new/raw-unormalised-data.csv")

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' Plot the stats
#' There is a large variance in sequencing depth. Some samples have been resequenced, hence
#' they will be merged in the next step
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6))
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' ## Combined samples
#' Combine them
count.table <- sapply(split.data.frame(t(count.table),sub("\\.[1,2]$","",colnames(count.table))),colSums)
count.stats <- sapply(split.data.frame(t(count.stats),sub("\\.[1,2]$","",colnames(count.stats))),colSums)
samples <- samples[!duplicated(sub("\\.[1,2]","",samples$SampleName)),c("SampleName","Line","LineTime","LineTimeSample","Time")]
samples$SampleName <- sub("\\.[1,2]","",samples$SampleName)
samples <- samples[match(colnames(count.table),samples$SampleName),]

#' Look at the stats
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' And plot
#' After merging the replicates, almost all samples have equal depth but OP42-2 - which presumably
#' failed the re-sequencing. T89-6 never made through to sequencing.
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6))
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is 79%
#' mean(unlist(count.stats["aligned",]/colSums(count.stats)))
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

#' Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The average cumulative coverage is around 100X
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' The observed distribution is relatively similar, a few samples show slight
#' diffences and there seem to be a trend of more spurious expression for
#' the T89 compared to the OP42. Probably an artefact of the alignment againt P. trichocarpa
plot.multidensity(lapply(1:ncol(count.table),function(i){log10(count.table[,i])}),
                  col=pal[as.integer(samples$Line)],
                  legend.x="topright",
                  legend=levels(samples$Line),
                  legend.col=pal[1:2],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = samples,
  design = ~ Line * Time)

dir.create(file.path("analysis","DE-new"),recursive = TRUE, showWarnings = FALSE)
save(dds,file="analysis/DE-new/DESeqDataset-20170703.rda")

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
write.csv(vst,"analysis/HTSeq-new/library-size-normalized_variance-stabilized_data.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd relationship. It is fairly linear,
#' meaning we can assume homoscedasticity. The slight initial trend / bump is
#' due to genes having few counts in a few subset of the samples and hence 
#' having a higher variability. This is expected.
meanSdPlot(vst[rowSums(count.table)>0,])

#' # QC on the normalised data
#' 
#' ## PCA
#' 
#' First perform a Principal Component Analysis (PCA) of the data
#'  to do a quick quality assessment; i.e. replicate should cluster
#' and the first 2-3 dimensions should be explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Plot the PCA 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Line)],
              pch=c(17,19)[as.integer(samples$Time)])
legend("topleft",pch=19,
       col=pal[1:2],
       legend=levels(factor(samples$Line)))

legend("topright",pch=c(17,19),
       legend=levels(factor(samples$Time)))
par(mar=mar)

#' Then the first two dimensions
#' The first two dimensions explains the time and the line
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Line)],
     pch=c(17,19)[as.integer(samples$Time)],
     main="Principal Component Analysis",sub="variance stabilized counts")

#' And the 2nd and 3rd dims
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples$Line)],
     pch=c(17,19)[as.integer(samples$Time)],
     main="Principal Component Analysis",sub="variance stabilized counts")

#' ## Heatmap
#' Select the genes - function
"geneSelect" <- function(cnt,splt,exp=1,nrep=2){
  rowSums(sapply(lapply(
    split.data.frame(t(cnt >= exp),splt)
    ,colSums), ">=", nrep)) >= 1
}

#' Define a cutoff
plot.multidensity(lapply(1:ncol(vst),function(i){vst[,i]}),
                  col=pal[as.integer(samples$LineTime)],
                  legend.x="topright",
                  legend=levels(samples$LineTime),
                  legend.col=pal[1:4],
                  legend.lwd=2,
                  main="sample vst expression distribution",
                  xlab="per gene vst expression (log2)")
abline(v=1:5,lty=2,col="grey")

#' We take a cutoff of 3 on the vst table. To keep a gene, it has to be expressed at least above a
#' vst value of 3 in at least 2 replicates of any condition (4 conditions total)
#' This is done only for plotting the heatmap. It could be used for filtering low spuriously expressed genes
#' for other analyses, but is not required for DE.
sel <- geneSelect(vst,samples$LineTime,exp=3)
message(sprintf("There are %s genes that pass the low spurious expression filter.", sum(sel)))

#' The vst clusters the sample first by Time and then by Line
#' The heatmap on the full data is squewed towards blue because of a few
#' highly expressed genes.
heatmap.2(vst[sel,],trace="none",labRow=FALSE,col=hpal)

#' Plot a dendrogram first
hc <- hclust(dist(t(vst[sel,])))
plot(hc)
plot(as.phylo.dendrogram(as.dendrogram(hc)))
plot(as.phylo.dendrogram(as.dendrogram(hc)),type="fan")

#' Saturate the expression - visualisation only
#' Select cutoffs
plot.multidensity(lapply(1:ncol(vst),function(i){vst[,i]}),
                  col=pal[as.integer(samples$LineTime)],
                  legend.x="topright",
                  legend=levels(samples$LineTime),
                  legend.col=pal[1:4],
                  legend.lwd=2,
                  main="sample vst expression distribution",
                  xlab="per gene vst expression (log2)")
abline(v=c(3,9),lty=2,col="grey")
vst.sat <- vst[sel,]
vst.sat[vst.sat < 3] <- 3
vst.sat[vst.sat > 9] <- 9
heatmap.2(vst.sat,trace="none",labRow=FALSE,col=hpal)

#' # Session Info
#' ``` {r session info, echo=FALSE}
#' sessionInfo()
#' ```
