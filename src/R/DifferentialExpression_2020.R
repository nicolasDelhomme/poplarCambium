#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Iryna Shutava"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/plotMA.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
#' ```{r edit1, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to change the variables in the 
#' plot to display the expression values accross your samples
#' The example below has 2 variables MGenotype and MDay. These 
#' need replacing by the variable(s) of interest in your project
#' ```
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=MDay,y=value,col=MGenotype,group=MGenotype)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,...){
    
    if(length(contrast)==1){
        res <- results(dds,name=contrast)
    } else {
        res <- results(dds,contrast=contrast)
    }
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
            
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
        
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
        message(sprintf("There are %s genes that are DE with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s",
                        sum(sel),
                        lfc,padj,expression_cutoff))
    }
    
    val <- rowSums(vst[sel,sample_sel])==0
    if (sum(val) >0){
        warning(sprintf("There are %s DE genes that have no vst expression in the selected samples",sum(val)))
        sel[sel][val] <- FALSE
    }    
    
    if(export){
        if(!dir.exists(default_dir)){
            dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
        }
        write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
        write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
    }
    if(plot){
        heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                  distfun = pearson.dist,
                  hclustfun = function(X){hclust(X,method="ward.D2")},
                  trace="none",col=hpal,labRow = FALSE,
                  labCol=labels[sample_sel],...
        )
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis/DE-new/DESeqDataset-20170703.rda"))
dds$LineTimeSample <- sub("\\.1$","",dds$LineTimeSample)

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE-new"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE-new/vst-aware.rda"))
write_csv(vst %>% as.data.frame %>% rownames_to_column("ID"),
          path=here("data/analysis/DE-new/vst-aware.csv"))

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' ## Results
#' ### T89 T1 vs. T0
contrast_T89 <- extract_results(dds=dds,vst=vst,contrast=c(0,0,1,1),
                              expression_cutoff=10,
                              default_prefix="T89-Time_T1_vs_T0",
                              labels=dds$LineTimeSample,
                              sample_sel=dds$Time %in% c("T0","T1") & dds$Line == "T89",
                              cexCol=0.8) 

#' ### OP42 T1 vs. T0
contrast_OP42 <- extract_results(dds=dds,vst=vst,contrast="Time_T1_vs_T0",
                              default_prefix="OP42-Time_T1_vs_T0",
                              labels=dds$LineTimeSample,
                              sample_sel=dds$Time %in% c("T0","T1") & dds$Line == "OP42",
                              cexCol=0.8) 

#' ### T89 vs. OP42 at T0
contrast_T0 <- extract_results(dds=dds,vst=vst,contrast="Line_T89_vs_OP42",
                               default_prefix="T0-T89_vs_OP42",
                               labels=dds$LineTimeSample, 
                               sample_sel=dds$Time == "T0",
                               cexCol=0.8) 

#' ### T89 vs. OP42 at T1
contrast_T1 <- extract_results(dds=dds,vst=vst,contrast=c(0,1,0,1),
                               default_prefix="T1-T89_vs_OP42",
                               labels=dds$LineTimeSample,
                               sample_sel=dds$Time == "T1",
                               cexCol=0.8)

#' ### Venn Diagram
#' Abbreviations:
#' 
#' * OP42: comparison T1 vs. T0 for OP42
#' 
#' * T89: comparison T1 vs. T0 for OP42
#' 
#' * T0: comparison T89 vs. OP42 at T0
#' 
#' * T1: comparison T89 vs. OP42 at T1
#' 
res.list <- list(OP42=contrast_OP42,
                 T1=contrast_T1,
                 T89=contrast_T89,
                 T0=contrast_T0)


#' #### All DE genes
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","all"),
                       NULL,
                       fill=pal[1:4]))

#' #### DE genes (up)
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","up"),
                      NULL,
                      fill=pal[1:4]))

#' #### DE genes (down)
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","dn"),
                       NULL,
                       fill=pal[1:4]))

#' ### Barplot
barplot(sapply(lapply(res.list,"[",-1),function(r){sapply(r,length)}) * c(1,-1),beside=TRUE)

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


