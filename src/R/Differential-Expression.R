#' ---
#' title: "Differential Expression in OP42 and T89 Cambium"
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
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(Glimma))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(VennDiagram))

#' Helper files
source("~/Git/UPSCb/src/R/plotMA.R")
source("~/Git/UPSCb/src/R/volcanoPlot.R")

#' Load saved data
load("analysis/DE-new/DESeqDataset-20170703.rda")

#' Setup graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' Cutoffs
fdr <- 0.01
lfc <- 0.5

#' Annotation
load("analysis/DE-new/annot.rda")

#' # Process
#' ## VST
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
write.csv(vst,"analysis/HTSeq-new/library-size-normalized_variance-stabilized-model-aware_data.csv")

#' ## DE 
#' Write a function
"doDE" <- function(dds){
  
  # Model and compute statistics
  dds <- DESeq(dds)
  
  # Plot the dispersion estimation
  plotDispEsts(dds)
  
  # Get the results
  res <- results(dds)
  
  # number of DE genes at a 1% FDR and min 0.5 LFC (Log2 Fold Change)
  sel <- res$padj <= fdr & !is.na(res$padj) & abs(res$log2FoldChange) >= lfc
  tab <- table(sign(res[sel,"log2FoldChange"]))
  message(sprintf("There are %s genes that are DE, %s down, %s up",sum(sel),tab[1],tab[2]))
  
  # plot the MA and volcano plot
  volcanoPlot(res,alpha=fdr)
  
  # Return object
  return(list(res,res[sel,]))
}

#' ## Differential Expression
#' ### OP42 Time
ddsOP42 <- dds[,grepl("OP",colnames(dds))]
design(ddsOP42) = ~Time
resOP42 <- doDE(ddsOP42)

#' Write the results
write.csv(resOP42[[1]],paste0("analysis/DE-new/OP42-T1-vs-T0-differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resOP42[[2]],
          paste0("analysis/DE-new/OP42-T1-vs-T0-differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))


#' ###  T89 Time
ddsT89 <- dds[,grepl("T89",colnames(dds))]
design(ddsT89) = ~Time

#'  We observe less DE genes in T89 compared to OP42. This is probably
#'  caused by the lack of one T0 replicate for T89. This makes DESeq2
#'  more conservative as can be observed by comparing the volcanoplot of
#'  T89 and OP42.
resT89 <- doDE(ddsT89)

#' Write the results
write.csv(resT89[[1]],paste0("analysis/DE-new/T89-T1-vs-T0-differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resT89[[2]],
          paste0("analysis/DE-new/T89-T1-vs-T0-differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))


#' ### Compare the results 
#' Plot the Venn diagram
#' 
#' The overlap is large (3/4 of the T89 genes are found in OP42), 
#' but as mentioned earlier the T89 DE analysis is presumbaly more 
#' constrained than the OP42 one (or in other words, less sensitive), so it is hard
#' to compare the results. To alleviate this issue, we need a more elaborate model
DE.list <- list(OP42=rownames(resOP42[[2]]),
                  T89=rownames(resT89[[2]]))

plot.new()
grid.draw(venn.diagram(DE.list,
                       filename=NULL,
                       col=pal[1:length(DE.list)],
                       category.names=names(DE.list)))

#' ## DE Line (Species) + Time
#' This means modeling the species effect and blocking it (i.e. ignoring it)
#' ### Model and compute statistics
design(dds) <- ~Line + Time

#'  We observe an almost similar number of gene differentially expressed as in OP42
#'  when we block the line factor. This is a good indication that this approach would
#'  make sense to compare T1 vs T0 in a species-agnostic manner. It seems to indicate
#'  that the processes are conserved in both species.
resLpT <- doDE(dds)

#' Write the results
write.csv(resLpT[[1]],paste0("analysis/DE-new/Line-plus-Time-T1-vs-T0-differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resLpT[[2]],paste0("analysis/DE-new/Line-plus-Time-T1-vs-T0-differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))

#' Create a Venn Diagram
DE.list <- list(OP42=rownames(resOP42[[2]]),
                T89=rownames(resT89[[2]]),
                LpT=rownames(resLpT[[2]]))
plot.new()
grid.draw(venn.diagram(DE.list, 
                       filename=NULL,
                       col=pal[1:length(DE.list)],
                       category.names=names(DE.list)))

#' ## DE Species * Time
#' Almost the same as above, but instead of "discarding" the Line effect, we
#' keep it and also look at the interaction between Line and Time. This assumes
#' that the variables Line and Time are not independent. 
design(dds) <- ~Line*Time
ddsLbT <- DESeq(dds) 

#' Plot the dispersion estimation
plotDispEsts(ddsLbT)

#' Check what results are available
resultsNames(ddsLbT)

#' ### Time effect different across Line
#'  We observe an almost similar number of gene differentially expressed as in OP42
#'  when we block the line faactor. This is a good indication that this approach would
#'  make sense to compare T1 vs T0 in a specie-agnostic manner. It seems to indicate
#'  that the processes are conserved in both species.
resLbT <- results(ddsLbT,name = "LineT89.TimeT1")
sel <- resLbT$padj <= fdr & !is.na(resLbT$padj) & abs(resLbT$log2FoldChange) > lfc
tab <- table(sign(resLbT[sel,"log2FoldChange"]))
message(sprintf("There are %s genes that are DE, %s down, %s up",sum(sel),tab[1],tab[2]))

#' Volcano plot
volcanoPlot(resLbT,alpha=0.01)

#' Heatmap
stopifnot(all(rownames(resLbT) == rownames(vst)))
heatmap.2(as.matrix(vst[sel,]),trace="none",labRow=FALSE,col=hpal)

#' ### Interactive plot
res.df <- as.data.frame(resLbT)
res.df$log10MeanNormCount <- log10(res.df$baseMean + 1)
idx <- rowSums(counts(ddsLbT)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
sel <- !is.na(resLbT$padj) & resLbT$padj <= fdr & abs(resLbT$log2FoldChange) >= lfc

#' When exploring this report, it is essential to remember that the 
#' dot plot on the left contain values that have simply been corrected
#' for their relative sequencing depth
glMDPlot(res.df,
         xval = "log10MeanNormCount",
         yval = "log2FoldChange",
         counts = counts(ddsLbT)[idx,],
         anno = annot[idx,],
         groups = ddsLbT$LineTime,
         samples = ddsLbT$SampleName,
         status = as.integer(sel[idx]),
         display.columns = c("GeneID","synonyms","Description","genelist_atg_id"),
         id.column = "GeneID",
         path = "analysis/DE-new",
         folder = "report-Line-Time-Interaction-Term-differential-expression",
         launch=FALSE)

#' Write the results
write.csv(resLbT,paste0("analysis/DE-new/Line-Time-Interaction-Term-differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resLbT[sel,],paste0("analysis/DE-new/Line-Time-Interaction-Term-differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))

#' ### A look at the Line
#' I.e. looking at genes that are different between lines time independent
resLine <- results(ddsLbT,name = "Line_T89_vs_OP42")
sel <- resLine$padj <= fdr & !is.na(resLine$padj) & abs(resLine$log2FoldChange) > lfc
tab <- table(sign(resLine[sel,"log2FoldChange"]))
message(sprintf("There are %s genes that are DE, %s down, %s up",sum(sel),tab[1],tab[2]))

#' Volcano plot
volcanoPlot(resLine,alpha=0.01)

#' Heatmap
heatmap.2(as.matrix(vst[sel,]),trace="none",labRow=FALSE,col=hpal)

#' ### Interactive plot
res.df <- as.data.frame(resLine)
res.df$log10MeanNormCount <- log10(res.df$baseMean + 1)
idx <- rowSums(counts(ddsLbT)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
sel <- !is.na(resLine$padj) & resLine$padj <= fdr & abs(resLine$log2FoldChange) >= lfc

#' When exploring this report, it is essential to remember that the 
#' dot plot on the left contain values that have simply been corrected
#' for their relative sequencing depth
glMDPlot(res.df,
         xval = "log10MeanNormCount",
         yval = "log2FoldChange",
         counts = counts(ddsLbT)[idx,],
         anno = annot[idx,],
         groups = ddsLbT$LineTime,
         samples = ddsLbT$SampleName,
         status = as.integer(sel[idx]),
         display.columns = c("GeneID","synonyms","Description","genelist_atg_id"),
         id.column = "GeneID",
         path = "analysis/DE-new",
         folder = "report-Line-Time_Line-effect_differential-expression",
         launch=FALSE)

#' Write the results
write.csv(resLine,paste0("analysis/DE-new/Line-Time_Line-effect_differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resLine[sel,],paste0("analysis/DE-new/Line-Time_Line-effect_differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))

#' ### A look at the Time
#' I.e. looking at genes that are different between lines time independent
resTime <- results(ddsLbT,name = "Time_T1_vs_T0")
sel <- resTime$padj <= fdr & !is.na(resTime$padj) & abs(resTime$log2FoldChange) > lfc
tab <- table(sign(resTime[sel,"log2FoldChange"]))
message(sprintf("There are %s genes that are DE, %s down, %s up",sum(sel),tab[1],tab[2]))

#' Volcano plot
volcanoPlot(resTime,alpha=0.01)

#' Heatmap
heatmap.2(as.matrix(vst[sel,]),trace="none",labRow=FALSE,col=hpal)

#' ### Interactive plot
res.df <- as.data.frame(resTime)
res.df$log10MeanNormCount <- log10(res.df$baseMean + 1)
idx <- rowSums(counts(ddsLbT)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
sel <- !is.na(resTime$padj) & resTime$padj <= fdr & abs(resTime$log2FoldChange) >= lfc

#' When exploring this report, it is essential to remember that the 
#' dot plot on the left contain values that have simply been corrected
#' for their relative sequencing depth
glMDPlot(res.df,
         xval = "log10MeanNormCount",
         yval = "log2FoldChange",
         counts = counts(ddsLbT)[idx,],
         anno = annot[idx,],
         groups = ddsLbT$LineTime,
         samples = ddsLbT$SampleName,
         status = as.integer(sel[idx]),
         display.columns = c("GeneID","synonyms","Description","genelist_atg_id"),
         id.column = "GeneID",
         path = "analysis/DE-new",
         folder = "report-Line-Time_Time-effect_differential-expression",
         launch=FALSE)

#' Write the results
write.csv(resTime,paste0("analysis/DE-new/Line-Time_Time-effect_differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resTime[sel,],paste0("analysis/DE-new/Line-Time_Time-effect_differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))

#' ## Only Line
#' ### Setup T0
ddsT0 <- dds[,grepl("T0",colData(dds)$Time)]
design(ddsT0) = ~Line
resT0 <- doDE(ddsT0)

#' Write the results
write.csv(resT0[[1]],paste0("analysis/DE-new/T0-T89-vs-OP42-differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resT0[[2]],paste0("analysis/DE-new/T0-T89-vs-OP42-differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))

#' ### Setup T1
ddsT1 <- dds[,grepl("T1",colData(dds)$Time)]
design(ddsT1) = ~Line
resT1 <- doDE(ddsT1)

#' Write the results
write.csv(resT1[[1]],paste0("analysis/DE-new/T1-T89-vs-OP42-differential-expression_FDR",fdr,"-LFC",lfc,".csv"))
write.csv(resT1[[2]],paste0("analysis/DE-new/T1-T89-vs-OP42-differential-expression_FDR",fdr,"-LFC",lfc,"-DE-only.csv"))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
