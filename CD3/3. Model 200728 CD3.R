require(treemap)
require(treemapify)
require(ggplot2)
require(plyr)
require(tidyverse)
require(dplyr)
require(reshape2)

# combining the data frame --------------------------------------------------------------

tol_CD3_lrt_trended = 0
mpe.a_CD3_lrt_trended = 0
mpe.r_CD3_lrt_trended = -1
sjs_CD3_lrt_trended = 1

lrt_trended <- glmLRT(fit.all.glmFit_trended,contrast=c(tol_CD3_lrt_trended,
                                       
                                        mpe.r_CD3_lrt_trended,
                                        mpe.a_CD3_lrt_trended,
                                        sjs_CD3_lrt_trended)) 

# top 10 tags --------------------------------------------------------------
topTags(lrt_trended, n=10) 

x.all <- as.data.frame(topTags(lrt_trended, adjust.method = "BH", n=dim(lrt_trended$table$PValue)))
x.all$Gene.ID <- rownames(x.all)

# plots --------------------------------------------------------------
edgeR_result <- as.data.frame(topTags(lrt_trended, adjust.method = "BH", n=dim(lrt_trended$table$PValue)))
dim(x.all)
topTags(lrt_trended, adjust.method = "BH", sort.by = "PValue")
summary(decideTests(lrt_trended, p=0.1, lfc=0.58))

pdf("pdf/MD_plot_actMPE.pdf", height=7, width=10)
plotMD(lrt_trended, cex = 0.5, hl.col = c("Green","red","black"),p=0.1)
abline(h=c(-0.58,0.58), col="blue", lwd = 2)
dev.off()

# Creating a df of the qlf DE genes --------------------------------------------------------------
x.all <- as.data.frame(topTags(lrt_trended, adjust.method = "BH", n=dim(lrt_trended$table$PValue)))
x.all_trended <- x.all
write.csv(x.all_trended[3],"csv/CD3.genelist.csv")
x.all_trended <- x.all_trended[c(3,9,5)]

head(x.all_trended)
names(x.all_trended) <- c("ID","Pvalue","logFC")
write.csv(x.all_trended,"csv/RSJS.RMPE.CD3.trended.csv",row.names = F)
