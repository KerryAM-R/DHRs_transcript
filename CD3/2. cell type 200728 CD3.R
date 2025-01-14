require(edgeR)
require(ggplot2)
library(reshape2)
library("gridExtra")
require(Biobase)
require(tidyverse)
library(magrittr)
library(tidyr)
require(ggrepel)

# dev.off()
names(mpe.2)
# combining the data frame #####
names(AMPE)[grep("CD3",names(AMPE))]
dim(tol.2[grep("CD3",names(tol.2))])
dim(sjs.2[grep("CD3",names(sjs.2))])
dim(RMPE[grep("CD3",names(RMPE))])
names(RMPE[grep("CD3",names(RMPE))])
all.a <- cbind(tol.2,RMPE,AMPE,sjs.2) 
names(all.a)
#all.a <- cbind(tol.2,sjs.2) 
all.a <- cbind(tol.2,mpe.2,sjs) 
dim(all.a)
all <-all.a[grep("CD3", names(all.a))]
names(all)
all_no_sex <- all
all_no_sex[(dim(all_no_sex)[2]+1)] <- rownames(all_no_sex)
names(all_no_sex)[(dim(all_no_sex)[2])] <- "Gene.ID"
all_no_sex_1 <- merge(d.anno,all_no_sex,by= "Gene.ID")
all_no_sex_2 <- all_no_sex_1[!(all_no_sex_1$Chrom=="X" | all_no_sex_1$Chrom=="Y"),]

row.names(all_no_sex_2) <- all_no_sex_2$Gene.ID

all_no_sex_2[all_no_sex_2$Gene.Name %in% "IL10",]
write.csv(all_no_sex_2[all_no_sex_2$Gene.Name %in% "IL10",],"IL10.CD3.csv")

write.table(all_no_sex_2,"csv/CD3_CD3.enrichment_no.sex.txt",row.names = F)
d.anno.sex <- all_no_sex_2[c(1:4)]
all <- all_no_sex_2[-c(1:4)]
head(all)
# creating the design matrix based on phenotype data #####

names.all.2 <- as.data.frame(do.call(rbind, strsplit(as.character(names(all)), "_")))
names.all.2
all_Cells <- factor(names.all.2$V4, levels = c("PBMC","CD3"))
all_Cells
patient <- factor(names.all.2$V2)
patient
batch <- factor(names.all.2$V3)
names.all.2[dim(names.all.2)[2]+1] <- paste(names.all.2$V3,names.all.2$V4, sep = "_")
batch_cell <- factor(names.all.2$V5)
batch_cell
names.all.2[dim(names.all.2)[2]+1] <- paste(names.all.2$V1,names.all.2$V4, sep = "_")
rash <- factor(names.all.2$V6, levels = c("Tolerant_PBMC","Tol_CD3","MPE.a_PBMC","MPE.a_CD3","MPE.r_PBMC","MPE.r_CD3","SJS_PBMC","RSJS_CD3"))
rash_cor <- factor(names.all.2$V1, levels = c("Tolerant","MPE","MPE.a","RSJS"))
names.all.2[dim(names.all.2)[2]+1] <- paste(names.all.2$V1,names.all.2$V3, sep = "_")
# DGEList object adn removing low expression based on counts per million cut-off
y.all <- DGEList(counts=all, genes = d.anno.sex)
AveLogCPM <- aveLogCPM(y.all)
nonzero <- rowSums(y.all$counts) >0
summary(nonzero)
y.all %<>% .[nonzero,]
y.all %<>% calcNormFactors
AveLogCPM <- aveLogCPM(y.all)
hist(AveLogCPM)
dev.off()
L <- mean(y.all$samples$lib.size) * 1e-6
M <- median(y.all$samples$lib.size) * 1e-6
c(L, M)
lcpm.cutoff <- log2(10/M + 2/L)
10/M + 2/L
lcpm.cutoff

lcpm_before <- cpm(y.all, log=TRUE)
dim(y.all$counts)
table(rowSums(y.all$counts==0)==28)

rash <- factor(names.all.2$V6, levels = c("Tol_PBMC","Tol_CD3","AMPE_PBMC","AMPE_CD3","MPE.r_PBMC","RMPE_CD3","RSJS_PBMC","RSJS_CD3"))
rash
lcpm.cutoff <- log2(10/M + 2/L) # determine counts per million cut-off
lcpm.cutoff
keep.exprs <- aveLogCPM(y.all) > lcpm.cutoff
table(keep.exprs)
2^lcpm.cutoff
y.all <- y.all[keep.exprs,, keep.lib.sizes=FALSE]
dim(y.all$counts)
lcpm_after <- cpm(y.all, log=TRUE)
unique(names.all.2$V6)
names.all.2$col <- with(names.all.2, ifelse(V6=="Tol_CD3","blue",
                                            ifelse(V6=="MPE.a_CD3","lightpink",
                                                   ifelse(V6=="MPE.r_CD3","darkorchid4","darkseagreen4"))))
# removing low counts #### -----
lcpm_melted <-melt(lcpm_before)
lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(lcpm_melted$Var2), "-")))
lcpm_melted$col.names <- lcpm_melted_names$V1

before <- ggplot(lcpm_melted, aes(x=value, colour=Var2))+
  geom_density()+ 
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 18)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=lcpm.cutoff, linetype="dashed", color = "grey") +
  theme(legend.position="bottom", legend.justification = "top",legend.title = element_blank())+
  theme(text=element_text(size=20,family="serif"),
        axis.title = element_text(colour="black", size=5,family="serif"),
        axis.text.x = element_text(colour="black",size=12,angle=90,hjust=0.6,vjust=.7,face="plain",family="serif"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
        legend.text=element_text(size=6),
        legend.position = "none")+
  guides(size=FALSE)+
  ggtitle("before low count filtering")

lcpm_melted <-melt(lcpm_after)
lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(lcpm_melted$Var2), "-")))

lcpm_melted$col.names <- lcpm_melted_names$V1

after <- ggplot(lcpm_melted, aes(x=value, colour=Var2))+
  geom_density()+ 
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 18)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=lcpm.cutoff, linetype="dashed", color = "grey") +
  theme(legend.position="bottom", legend.justification = "top",legend.title = element_blank())+
  theme(text=element_text(size=20,family="serif"),
        axis.title = element_text(colour="black", size=5,family="serif"),
        axis.text.x = element_text(colour="black",size=12,angle=90,hjust=0.6,vjust=.7,face="plain",family="serif"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
        legend.text=element_text(size=6),
        legend.position = "none")+
  guides(size=FALSE)+
  ggtitle("After low count filtering")

pdf("pdf/count.filtering.pdf",height=6,width=12)
grid.arrange(before,after,ncol=2) # printing out the graphs
dev.off()

# normalisation -----
rm(before)
rm(after)
y.all$samples$norm.factors <- 1
y.all$samples$norm.factors
lcpm_un_norm <- cpm(y.all, log=TRUE)

lcpm_melted <-melt(lcpm_un_norm)
lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(lcpm_melted$Var2), "_")))
lcpm_melted$col.names <- paste(lcpm_melted_names$V1,lcpm_melted_names$V5,sep="_")

before <- ggplot(lcpm_melted, aes(x=Var2,y=value, fill=col.names))+
  geom_boxplot()+ 
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 18)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  theme(legend.position="bottom", legend.justification = "top",legend.title = element_blank())+
  theme(text=element_text(size=20,family="serif"),
        axis.title = element_text(colour="black", size=5,family="serif"),
        axis.text.x = element_text(colour="black",size=12,angle=90,hjust=0.6,vjust=.7,face="plain",family="serif"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
        legend.text=element_text(size=12),
        legend.position = "none")+
  guides(size=FALSE)+
  scale_x_discrete(labels=names.all.2$V8)


y.all <- calcNormFactors(y.all, method = "TMM")
y.all$samples$norm.factors
lcpm_norm <- cpm(y.all, log=TRUE)

lcpm_melted <-melt(lcpm_norm)
lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(lcpm_melted$Var2), "_")))
lcpm_melted$col.names <- paste(lcpm_melted_names$V1,lcpm_melted_names$V5,sep="_")

after <- ggplot(lcpm_melted, aes(x=Var2,y=value, fill=col.names))+
  geom_boxplot()+ 
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  theme_bw(base_size = 18)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") +
  theme(legend.position="bottom", legend.justification = "top",legend.title = element_blank())+
  theme(text=element_text(size=20,family="serif"),
        axis.title = element_text(colour="black", size=5,family="serif"),
        axis.text.x = element_text(colour="black",size=12,angle=90,hjust=0.6,vjust=.7,face="plain",family="serif"),
        axis.text.y = element_text(colour="black",size=16,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(colour="black",size=16,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
        legend.text=element_text(size=12),
        legend.position = "none")+
  guides(size=FALSE)+
  scale_x_discrete(labels=names.all.2$V8)

pdf("pdf/normalisation.pdf",heigh=6,width =10)
grid.arrange(before,after,ncol=2)
dev.off()

# PCA analysis for statiscial outliers -----
log <- cpm(y.all, log=T)
head(log)
project.pca <- prcomp(t(log))
head(project.pca$x)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

pdf("pdf/Scree_plot.pdf",heigh=10,width =10)
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
dev.off()


round(project.pca.proportionvariances, 2)

head(names.all.2)
All_col <- factor(names.all.2$V1)
pairs(project.pca$x[,1:5], col=as.numeric(All_col), main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
pairs(project.pca$x[,6:10], col=as.numeric(All_col), main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()
#Plots scatter plot for PC 1 and 2
a <- as.data.frame(project.pca$x)
head(names.all.2)
a[dim(a)[2]+1] <- names.all.2$V1
a[dim(a)[2]+1] <- names.all.2$V2
a[dim(a)[2]+1] <- names.all.2$V3
a[dim(a)[2]+1] <- names.all.2$V4
a[dim(a)[2]+1] <- names.all.2$V5
a[dim(a)[2]+1] <- names.all.2$V6
a[dim(a)[2]+1] <- names.all.2$V7
a <- as.data.frame(a)
(dim(a)[2]-7):dim(a)[2]
names(a)[(dim(a)[2]-6):dim(a)[2]] <- c("rash","ID","batch","cells","batch_cells","rash_cells","rash_batch")

require(ggplot2)
#dev.off()
require("ggrepel")
# a$ID_rash[order(a$ID_rash)]
library(extrafont)
loadfonts(device = "pdf")

pdf("pdf/PCA_1_2.pdf", height=5, width=6)

ggplot(a, aes(x = PC1, y = PC2, colour = rash, label = ID)) +
  geom_point(size=4) +
  theme_bw() +
  scale_x_continuous(
    name = paste("PC1,", round(project.pca.proportionvariances[1], 2), "%"), 
    # limits = c(-80, 120)
  ) +
  scale_y_continuous(
    name = paste("PC2,", round(project.pca.proportionvariances[2], 2), "%")
  ) +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.justification = "top",
    text = element_text(family = "Times New Roman"), # Apply Times New Roman
    axis.title = element_text(family = "Times New Roman"), 
    axis.text = element_text(family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman"),
    legend.title = element_text(family = "Times New Roman")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values=c("#e37c1d","#0076c0","#a30234","grey")) +
  geom_text_repel(
    segment.alpha = 0.25, 
    show.legend = F,
    col = "black",
    box.padding = unit(0.5, 'lines'), 
    max.overlaps = Inf,
    family = "Times New Roman" # Apply font to geom_text_repel
  ) 

dev.off()

# mean(a$PC2) +3*sd(a$PC2)

# openPDF("pdf/PCA_1_2.pdf")
pdf("pdf/PCA_3-4.pdf", height=7, width=10)
ggplot(a, aes(x = PC3, y = PC4, colour = rash)) +
  geom_point(size=4) +
  theme_bw()  +
  scale_x_continuous(name = paste("PC3,", round(project.pca.proportionvariances[3], 2), "%"))+
  scale_y_continuous(name = paste("PC4,", round(project.pca.proportionvariances[4], 2), "%")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(legend.justification = "top") +
  scale_color_discrete("group") +
  #geom_text(aes(x = PC3+10, y = PC4,label=ID),size = 3.5, col = "black")+
  guides(colour = guide_legend(override.aes = list(size = 5))) 

dev.off()

######
# Design matrix
#####
names.all.2

rash <- factor(names.all.2$V6, levels = c("Tol_CD3","RMPE_CD3","AMPE_CD3","RSJS_CD3"))

y.all


design.all <- model.matrix(~0+rash)
design.all
y.all <- estimateDisp(y.all, design=design.all)

logFC <- predFC(y.all,design.all,prior.count=1)
pdf("pdf/BVC.pdf", height=7, width=10)
plotBCV(y.all)
dev.off()

fit.all.glmFit_trended <- glmQLFit(y.all,design.all, dispersion = y.all$trended.dispersion,robust=T)
plotQLDisp(fit.all.glmFit_trended)
write.table(y.all$genes[1],"csv/CD3.list.csv",row.names = F,col.names = F,sep=",",quote = F)
