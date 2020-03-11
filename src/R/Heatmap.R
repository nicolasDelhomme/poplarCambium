#' ---
#' title: "Gene of interest heatmap"
#' author: "Nicolas Delhomme & Alok Ranjan"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' # Environment
#' Set the working dir
setwd("/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium")
#' ```

#' Libs
suppressPackageStartupMessages(library(gplots))

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Gene lists
aux <- unique(read.delim("~/Git/UPSCb/projects/bellini_cambium/doc/Aux-gene-list.txt",as.is=TRUE))
ctk <- unique(read.delim("~/Git/UPSCb/projects/bellini_cambium/doc/Cytokinin-gene-list.txt",as.is=TRUE))
ja <- unique(read.delim("~/Git/UPSCb/projects/bellini_cambium/doc/JA-gene-list.txt",as.is=TRUE))
et <- unique(read.delim("~/Git/UPSCb/projects/bellini_cambium/doc/ET-gene-list.txt",as.is=TRUE))

#' Data
vst <- read.csv("analysis/HTSeq-new/library-size-normalized_variance-stabilized-model-aware_data.csv",row.names = 1)
rownames(vst) <- sub("\\.0$","",rownames(vst))

#' Sample info
samples <- read.csv("~/Git/UPSCb/projects/bellini_cambium/doc/samples.csv")
samples <- samples[!duplicated(sub("\\.[1,2]$","",samples$SampleName)),]
samples$SampleName <- sub("-",".",sub("\\.[1,2]$","",samples$SampleName))
samples <- samples[match(colnames(vst),samples$SampleName),]

#' Out dir
dir.create("analysis/Heatmaps",showWarnings = FALSE)

#' # Heatmap
#' ## Auxin
aux.dat <- as.matrix(vst[aux[,1],])
stopifnot(all(rownames(aux.dat) == aux))
write.csv(aux.dat,"analysis/Heatmaps/Auxin-genes_vst-expression-data.csv")
aux.mean <- sapply(split.data.frame(t(aux.dat),samples$LineTime),colMeans)
sel <- featureSelect(aux.mean,factor(colnames(aux.mean)),exp = 6,nrep = 1)

pdf("analysis/Heatmaps/Auxin-genes_vst-expression_heatmap.pdf",width=8,height=12)
heatmap.2(aux.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")
dev.off()

#' ## JA
ja.dat <- as.matrix(vst[ja[,1],])
stopifnot(all(rownames(ja.dat) == ja))
write.csv(ja.dat,"analysis/Heatmaps/ja-genes_vst-expression-data.csv")
ja.mean <- sapply(split.data.frame(t(ja.dat),samples$LineTime),colMeans)
sel <- featureSelect(ja.mean,factor(colnames(ja.mean)),exp = 1,nrep = 1)
heatmap.2(ja.mean[sel,],col=hpal,trace="none",margins=c(8,8),Colv=FALSE,dendrogram = "row")
sel <- featureSelect(ja.mean,factor(colnames(ja.mean)),exp = 4,nrep = 1)
heatmap.2(ja.mean[sel,],col=hpal,trace="none",margins=c(8,8),Colv=FALSE,dendrogram = "row")

pdf("analysis/Heatmaps/ja-genes_vst-expression_heatmap.pdf",width=8,height=12)
heatmap.2(ja.mean[sel,],col=hpal,trace="none",margins=c(8,8),Colv=FALSE,dendrogram = "row")
dev.off()

#' ## cytokinin
ctk.dat <- as.matrix(vst[ctk[,1],])
stopifnot(all(rownames(ctk.dat) == ctk))
write.csv(ctk.dat,"analysis/Heatmaps/ctk-genes_vst-expression-data.csv")
ctk.mean <- sapply(split.data.frame(t(ctk.dat),samples$LineTime),colMeans)
sel <- featureSelect(ctk.mean,factor(colnames(ctk.mean)),exp = 1,nrep = 1)
heatmap.2(ctk.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")
sel <- featureSelect(ctk.mean,factor(colnames(ctk.mean)),exp = 4,nrep = 1)
heatmap.2(ctk.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")

pdf("analysis/Heatmaps/ctk-genes_vst-expression_heatmap.pdf",width=8,height=12)
heatmap.2(ctk.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")
dev.off()

#' ## Ethylene
et.dat <- as.matrix(vst[et[,1],])
stopifnot(all(rownames(et.dat) == et))
write.csv(et.dat,"analysis/Heatmaps/et-genes_vst-expression-data.csv")
et.mean <- sapply(split.data.frame(t(et.dat),samples$LineTime),colMeans)
sel <- featureSelect(et.mean,factor(colnames(et.mean)),exp = 1,nrep = 1)
heatmap.2(et.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")
sel <- featureSelect(et.mean,factor(colnames(et.mean)),exp = 4,nrep = 1)
heatmap.2(et.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")

pdf("analysis/Heatmaps/et-genes_vst-expression_heatmap.pdf",width=8,height=12)
heatmap.2(et.mean[sel,],col=hpal,trace="none",margins=c(8,12),Colv=FALSE,dendrogram = "row")
dev.off()


# heatmap.2(aux.mean[sel,],col=hpal,trace="none",Colv = FALSE,dendrogram = "row",RowSideColors = c(rep("2",13),rep("3",20)))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
