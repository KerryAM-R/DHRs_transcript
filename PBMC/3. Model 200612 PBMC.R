require(treemap)
require(treemapify)
require(ggplot2)
require(plyr)
require(tidyverse)
require(dplyr)
require(reshape2)

# combining the data frame --------------------------------------------------------------

tol_PBMC_lrt_trend = 0
mpe.a_PBMC_lrt_trend = -1
mpe.r_PBMC_lrt_trend = 0
sjs_PBMC_lrt_trend = 1

# removed batch from model as it was skewing the results...
lrt_trend <- glmLRT(fit.all.glmFit_trended,contrast=c(tol_PBMC_lrt_trend,
                                        mpe.a_PBMC_lrt_trend,
                                        mpe.r_PBMC_lrt_trend,
                                        sjs_PBMC_lrt_trend)) 

 # top 10 tags --------------------------------------------------------------
topTags(lrt_trend, n=10) 

x.all <- as.data.frame(topTags(lrt_trend, adjust.method = "BH", n=dim(lrt_trend$table$PValue)))
#write.csv(x.all,"resMPE.all.csv")
x.all_trended <- as.data.frame(topTags(lrt_trend, adjust.method = "BH", n=dim(lrt_trend$table$PValue)))
names(x.all_trended)[c(5:9)] <- paste(names(x.all_trended)[c(5:9)],"Trended",sep="_")
x.all$Gene.ID <- rownames(x.all)
head(x.all)

# plots --------------------------------------------------------------
edgeR_result <- as.data.frame(topTags(lrt_trend, adjust.method = "BH", n=dim(lrt_trend$table$PValue)))

topTags(lrt_trend, adjust.method = "BH", sort.by = "PValue")
summary(decideTests(lrt_trend, p=0.1,lfc = 0.58))
lrt_trend$comparison <- gsub("rash","",lrt_trend$comparison)

pdf("pdf/plot_MD.SJS.rMPE.pdf", height=7, width=10)
plotMD(lrt_trend, cex = 0.5, hl.col = c("Green","red","black"),p=0.05)
abline(h=c(-0.58,0.58), col="blue", lwd = 2)
dev.off()

# Creating a df of the qlf DE genes --------------------------------------------------------------
x.all <- as.data.frame(topTags(lrt_trend, adjust.method = "BH", n=dim(lrt_trend$table$PValue)))
x.all$count <- 1

x.all_trended <-x.all
x.all_trended <- x.all_trended[c(3,9,5)]
head(x.all_trended)
names(x.all_trended) <- c("ID","Pvalue","logFC")
write.csv(x.all_trended,"csv/SJS.aMPE.PBMC.trended.csv",row.names = F)
