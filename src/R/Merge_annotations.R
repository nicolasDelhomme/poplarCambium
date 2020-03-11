library(parallel)
setwd("/mnt/picea/projects/aspseq/cbellini/24_Bellini_cambium/")

popgenie_annaotation <- read.csv("analysis/DE/gene_table.csv",as.is = TRUE)

annot <- unique(popgenie_annaotation[,-1])

annot <- do.call(rbind,mclapply(unique(annot$GeneNames),function(g,an){
  r<-an[an$GeneNames == g,]
  if(nrow(r)>1){
    r[1,4] <- paste(unique(r[,4]),collapse=";")
    r[1,5] <- paste(unique(r[,5]),collapse=";")
    r <- r[1,]
  }
  return(r)
},annot,mc.cores=16))

#T1-T89-vs-op42 annotation
results_csv<-read.csv("analysis/DE/T1-T89-vs-OP42-differential-expression.csv",row.names = 1)
res <- cbind(results_csv,annot[match(sub("\\.0$","",rownames(results_csv)),annot$GeneNames),])
write.csv(res,file="analysis/DE/T1-T89-vs-OP42-differential-expression_annotated.csv")

#T0-T89-vs-op42 annotation
results_csv<-read.csv("analysis/DE/T0-T89-vs-OP42-differential-expression.csv",row.names = 1)
popgenie_annaotation <- read.csv("analysis/DE/gene_table.csv",as.is = TRUE)
res2 <- cbind(results_csv,annot[match(sub("\\.0$","",rownames(results_csv)),annot$GeneNames),])
write.csv(res2,file="analysis/DE/T0-T89-vs-OP42-differential-expression_annotated.csv")

##T89-T1-vs-T0

results_csv<-read.csv("analysis/DE/T89-T1-vs-T0-differential-expression.csv",row.names = 1)
popgenie_annaotation <- read.csv("analysis/DE/gene_table.csv",as.is = TRUE)
res3 <- cbind(results_csv,annot[match(sub("\\.0$","",rownames(results_csv)),annot$GeneNames),])
write.csv(res3,file="analysis/DE/T89-T1-vs-T0-differential-expression_annotated.csv")

##op42-T1-vs-T0
results_csv<-read.csv("analysis/DE/OP42-T1-vs-T0-differential-expression.csv",row.names = 1)
popgenie_annaotation <- read.csv("analysis/DE/gene_table.csv",as.is = TRUE)
res4 <- cbind(results_csv,annot[match(sub("\\.0$","",rownames(results_csv)),annot$GeneNames),])
write.csv(res4,file="analysis/DE/OP42-T1-vs-T0-differential-expression_annotated.csv")


