#! /usr/bin/env Rscript

packrat::init("~/runs/eyu8/library/panelCNMOPS")
setwd("~/runs/eyu8/data/CNV/panelcnMOPS")

library(panelcn.mops)
library(tidyverse)
library(glue)
library(dplyr)

normType <- "quant"
sizeFactor <- "quant"
CN1 <- 0.57
CN2 <- 1
CN3 <- 1.46
priorImpact <- 1
version <- 1


duplicates <- readLines("../../MIP/txt/data/NeuroX_duplicate_v3.data")
duplicates <- paste0('MIP.',duplicates,'.clean.bam',sep='')


#poorSamples <- readLines("../MIP/txt/files/PARK2_bad.csv")
#poorSamples <- paste0('MIP.',poorSamples,'.clean.bam',sep='')

#good_gene <- paste(readLines("../MIP/txt/cov/good_gene_above_90.cov"),collapse="|")
good_gene <- readLines("../../MIP/txt/cov/good_gene_above_90.cov")

MLPA_cnv <- paste0('MIP.',readLines("txt/MLPA_cnv"),'.clean.bam',sep='')
MLPA_wild <- paste0('MIP.',readLines("txt/MLPA_wild"),'.clean.bam',sep='')


countWindows <- getWindows("named_mip_target_genes.bed")
bed <- countWindows[countWindows$gene %in% good_gene,]


pd_number <- c("pd_15","pd_19","pd_22","pd_26","pd_28","pd_55","pd_65","pd_71","pd_76")
ctrl_number <- c("control_19","control_22","control_26","control_28","control_55","control_65","control_71","control_76")



pd.list <- lapply(pd_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
#for(i in 1:length(pd.list)){names(ranges(pd.list[[i]])) <- rownames(countWindows)}
pd_PARK2 <-  lapply(pd.list, function(x) x[as.numeric(rownames(bed)),])
pd.tmp1 <- lapply(pd_PARK2, function(x) x[ , !(names(mcols(x)) %in% duplicates)])
pd.old.dafr <- lapply(pd.tmp1, function(x) x[ , !(names(mcols(x)) %in% duplicates)])
all_pd_2  <- Reduce(cbind,lapply(2:9,function(x) elementMetadata(pd.old.dafr[[x]])))
all_pd <- pd.old.dafr[[1]]
elementMetadata(all_pd) <- cbind(elementMetadata(all_pd),all_pd_2)

ctrl.list <- lapply(ctrl_number, function(x)readRDS(paste("raw/",x,"_all.rds",sep="")))
#for(i in 1:length(ctrl.list)){names(ranges(ctrl.list[[i]])) <- countWindows$gene}
#ctrl_PARK2 <-  lapply(ctrl.list, function(x) x[grep(names(x), pattern = good_gene),])
ctrl_PARK2 <-  lapply(ctrl.list, function(x) x[as.numeric(rownames(bed)),])
ctrl.tmp1 <- lapply(ctrl_PARK2, function(x) x[ , !(names(mcols(x)) %in% duplicates)])
ctrl.dafr <- lapply(ctrl.tmp1, function(x) x[ , !(names(mcols(x)) %in% duplicates)])
all_ctrl_2  <- Reduce(cbind,lapply(2:8,function(x) elementMetadata(ctrl.dafr[[x]])))
all_ctrl <- ctrl.dafr[[1]]
elementMetadata(all_ctrl) <- cbind(elementMetadata(all_ctrl),all_ctrl_2)

all_samples <- all_pd
elementMetadata(all_samples) <- cbind(elementMetadata(all_pd),elementMetadata(all_ctrl) )


CNV_samples <- all_samples[ , names(mcols(all_samples)) %in% MLPA_cnv]
noCNV_samples <- all_samples[ , names(mcols(all_samples)) %in% MLPA_wild]

ctrl.ref.list <- lapply(ctrl.dafr, function(x) x[ , !(names(mcols(x)) %in% c(MLPA_cnv,MLPA_wild))])

ctrl.ref <- Reduce(cbind,lapply(ctrl.ref.list,elementMetadata))


#loop through different test

getResult <- function(test,I,normType,sizeFactor,priorImpact){
    XandCB <- test
    elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), ctrl.ref)

    resultlist <- runPanelcnMops(XandCB, testiv = 1:ncol(elementMetadata(test)),countWindows = bed, I = as.numeric(I), normType = normType, sizeFactor = sizeFactor, priorImpact = as.numeric(priorImpact))
    sampleNames <- colnames(elementMetadata(test))
    return(createResultTable(resultlist = resultlist, XandCB = XandCB, countWindows = countWindows, sampleNames = sampleNames))
}

I <- c(0.025, CN1, CN2, CN3, 2)
resulttable <- getResult(CNV_samples,I,normType,sizeFactor,priorImpact)
score_list <- lapply(1:length(resulttable), function(n){
                         result <- resulttable[[n]]
                         PARK2 <- result[result$Gene == "PARK2",]
                         TP_list <- ifelse(FALSE %in% (PARK2$CN == "CN2"), 1, 0)
                         return(c(levels(PARK2$Sample),TP_list))
})

s <- as.data.frame(Reduce(rbind,score_list))
colnames(s) <- c("Sample","SCORE")

write.csv(s,"PD_cases_performance.csv",row.names=FALSE)

resulttable <- getResult(noCNV_samples,I,normType,sizeFactor,priorImpact)
score_list <- lapply(1:length(resulttable), function(n){
                         result <- resulttable[[n]]
                         PARK2 <- result[result$Gene == "PARK2",]
                         TN_list <- ifelse(FALSE %in% (PARK2$CN == "CN2"), 0, 1)
                         return(c(levels(PARK2$Sample),TN_list))
})
s <- as.data.frame(Reduce(rbind,score_list))
colnames(s) <- c("Sample","SCORE")
write.csv(s,"PD_controls_performance.csv",row.names=FALSE)

#F1 <- 2*s[1]/(2*s[1] + s[2] + score[2])
#message("Accuracy:",s[1]/(s[1] + s[2]))
#FPR <- score[2]/(score[1] + score[2]) 
#message("FPR:",FPR)
#message("Out of:",score[1] + score[2])
#cat(c(normType, sizeFactor, CN1, CN2, CN3, priorImpact, F1, FPR),file = glue('param/test_{version}'))
#message(glue("Wrote param/test_{version}"))
#message("TP:",s[1],"FP:",score[2],"TN:",score[1],"FN:",s[2])
#                        }               
#                    }
#                }
#            }
#        }
#    }



#end of loop

