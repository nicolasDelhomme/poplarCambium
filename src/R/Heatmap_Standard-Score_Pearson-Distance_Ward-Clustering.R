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
suppressPackageStartupMessages(library(hyperSpec))

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Gene lists
arf <- unique(read.delim("~/Git/UPSCb/projects/bellini_cambium/doc/NewHeatmapAnalysis/Auxin-ARF.txt",
                         as.is=TRUE,header=FALSE))

#' Data
vst <- read.csv("analysis/HTSeq-new/library-size-normalized_variance-stabilized-model-aware_data.csv",row.names = 1)
rownames(vst) <- sub("\\.0$","",rownames(vst))

#' Sample info
samples <- read.csv("~/Git/UPSCb/projects/bellini_cambium/doc/samples.csv")
samples <- samples[!duplicated(sub("\\.[1,2]$","",samples$SampleName)),]
samples$SampleName <- sub("-",".",sub("\\.[1,2]$","",samples$SampleName))
samples <- samples[match(colnames(vst),samples$SampleName),]

#' Out dir
dir.create("analysis/NewHeatmaps",showWarnings = FALSE)

#' Function
".heatmapMe" <- function(vst,goi,nam,exp=1){
  stopifnot(all(goi[,1] %in% rownames(vst)))
  dat <- as.matrix(vst[goi[,1],])
  stopifnot(all(rownames(dat) == goi))
  write.csv(dat,paste0("analysis/NewHeatmaps/",nam,"_genes_vst-expression_cutoff-expression-",exp,"_data.csv"))
  dat.mean <- sapply(split.data.frame(t(dat),samples$LineTime),colMeans)
  sel <- featureSelect(dat.mean,factor(colnames(dat.mean)),exp = exp,nrep = 1)
  
  # expression
  heatmap.2(dat[sel,],
            col=hpal,
            trace="none",margins=c(8,12),
            Colv=FALSE,dendrogram = "row",
            main=paste0(nam," expression"))
  
  # mean expression
  heatmap.2(dat.mean[sel,],
            col=hpal,
            trace="none",margins=c(8,12),
            Colv=FALSE,dendrogram = "row",
            main=paste0(nam," mean expression"))
  
  # z-score
  s.dat <- t(scale(t(dat[sel,])))
  heatmap.2(s.dat,distfun = pearson.dist,
            hclustfun = function(X){hclust(X,method="ward.D")},
            col=hpal,
            trace="none",margins=c(8,12),
            Colv=FALSE,dendrogram = "row",
            main=paste0(nam," z-score"))
  
  # mean z-score
  s.dat.mean <- t(scale(t(dat.mean[sel,])))
  heatmap.2(s.dat.mean,distfun = pearson.dist,
            hclustfun = function(X){hclust(X,method="ward.D")},
            col=hpal,
            trace="none",margins=c(8,12),
            Colv=FALSE,dendrogram = "row",
            main=paste0(nam," mean z-score"))
  
  #,
            #colsep=1:ncol(s.dat.mean),
            #  rowsep=1:nrow(s.dat.mean),
            #  sepcolor="black")
  
  
}

#' # Heatmap
#' ## Auxin ARFs
pdf("analysis/Heatmaps/Auxin-ARF-genes_vst-expression_heatmap.pdf",width=8,height=12)
.heatmapMe(vst,arf,"Auxin-ARF",exp=2)
.heatmapMe(vst,arf,"Auxin-ARF",exp=4)
.heatmapMe(vst,arf,"Auxin-ARF",exp=8)
dev.off()

#' ## Auxin IAA
pdf("analysis/Heatmaps/Auxin-ARF-genes_vst-expression_heatmap.pdf",width=8,height=12)
.heatmapMe(vst,,"Auxin-",exp=2)
.heatmapMe(vst,,"Auxin-",exp=4)
.heatmapMe(vst,,"Auxin-",exp=8)
dev.off()



#' ## ...


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
