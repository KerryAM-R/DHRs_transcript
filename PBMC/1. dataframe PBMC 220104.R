#BiocManager::install("edgeR")
library(edgeR)
library(dplyr)
library(pheatmap)
##### Master identification sheet -----
Master.ID <- read.csv("../master_ID_updated.ID.csv", header = T)
Master.ID

# head(mpe.2[grep("DR20031",names(mpe.2))])

#read in dataframe and replace names -----
df <- read.table("../NonStrandedCounts-withNames.txt",sep='\t',header=T,stringsAsFactors = F)

#df <- read.table("added TR2_3_5.txt",sep='\t',header=T,stringsAsFactors = F)
head(df)
names(df) <- gsub('X','',names(df))
rownames(df) <- df$Gene.ID
df.names <- as.data.frame(names(df))
df.names
names(df.names) <- "oldPrefix"
ID.df <- read.csv("../new_prefex.csv", header = T)
head(ID.df)

# ID.df$Cell <- gsub("CD3","CD3",ID.df$Cell)

df.names.all <- merge(df.names,ID.df, by ="oldPrefix", sort = F) 
df.names.all
df.names.all.master <- merge(df.names.all,Master.ID,by="ID", sort = F)
df.names.all.master.2 <- merge(df.names,df.names.all.master,by="oldPrefix")
df.names.all.master.2
df.names.all.master.2$replacement_ID <- paste(df.names.all.master.2$abbreviated.ID,df.names.all.master.2$no., df.names.all.master.2$Batch,df.names.all.master.2$Cell, sep ="_")

head(df.names.all.master.2)

write.csv(df.names.all.master.2,"csv/df.names.all.master.2.csv")
d.anno <- df[,c(1:4)]

#write.csv(d.anno,"d.anno.csv")

head(d.anno)
names(d.anno) <- gsub('X','',names(d.anno))
d.data <- df[,-c(1:4)] #63677
names(d.data)
names(d.data) <- df.names.all.master.2$replacement_ID
names(d.data)
head(d.data)
names(df)
#groups -----

sjs <-d.data[grep("SJS", names(d.data))]
tol <-d.data[grep("Tol", names(d.data))]

names(tol)[grep("_11_",names(tol))]

names(tol)
mpe <-d.data[grep("MPE", names(d.data))]
unique(names(mpe))
names(sjs)[grep("_15_",names(sjs))]

names(sjs[grep("PBMC", names(sjs))])
names(mpe[grep("PBMC", names(mpe))])
names(mpe[grep("CD3", names(mpe))])
names(sjs)

# remove tolerant outliers -----
names(tol)
names(tol[grep("CD3", names(tol))])
names(tol)
tol.1 <- tol
names(tol.1)
tol.2 <- tol.1[,!names(tol.1) %in% c("Tol_4_B2_CD3","Tol_1_B2_PBMC")]
names(tol.2)

names(sjs)
sjs.1 <- sjs
round((sjs$RSJS_4_B1_PBMC+sjs$RSJS_4_B3_PBMC)/2,0)

sjs.1$RSJS_4av_B1_PBMC <- round((sjs$RSJS_4_B1_PBMC+sjs$RSJS_4_B3_PBMC)/2,0)
names(sjs.1)
sjs.2<-  sjs.1[,!names(sjs.1) %in% c("RSJS_4_B1_PBMC","RSJS_4_B3_PBMC")]
dim(sjs.2)
names(sjs.2)
# remove mpe outlier/tech

mpe.1 <- mpe
names(mpe.1)

names(mpe[grep("PBMC", names(mpe))])
names(mpe[grep("CD3", names(mpe))])

mpe.1$RMPE_3av_B1_PBMC <- round((mpe.1$RMPE_3_B1_PBMC+mpe.1$RMPE_3_B1_PBMC.1)/2,0) # changed from 4 to 3 to match the functional data
mpe.1$RMPE_4av_B1_PBMC <- round((mpe.1$RMPE_4_B1_PBMC+mpe.1$RMPE_4_B1_PBMC.1)/2,0) # changed from 5 to 4 to match the functional data
mpe.1$RMPE_5av_B3_PBMC <- round((mpe.1$RMPE_5_B3_PBMC+mpe.1$RMPE_5_B1_PBMC)/2,0) # changed from 8 to 5 to match the functional data

mpe.1$RMPE_3av_B2_CD3 <- round((mpe.1$RMPE_3_B2_CD3+mpe.1$RMPE_3_B3_CD3)/2,0) # changed from 4 to 3 to match the functional data
mpe.1$RMPE_7av_B3_CD3 <- round((mpe.1$RMPE_7_B2_CD3+mpe.1$RMPE_7_B3_CD3)/2,0) # changed from 10 to 7 to match the functional data

mpe.2 <- mpe.1[,!names(mpe.1) %in% c("RMPE_3_B1_PBMC","RMPE_3_B1_PBMC.1",
                                     "RMPE_3_B2_CD3","RMPE_3_B3_CD3",
                                     "RMPE_4_B1_PBMC","RMPE_4_B1_PBMC.1",
                                     "AMPE_3_B3_PBMC",
                                     "RMPE_5_B3_PBMC","RMPE_5_B1_PBMC",
                                     "RMPE_7_B2_CD3","RMPE_7_B3_CD3"
                                     )]
# active and resolved mpe cases
RMPE <- mpe.2[grep("RMPE",names(mpe.2))]
names(RMPE)




AMPE <- mpe.2[grep("AMPE",names(mpe.2))]
names(AMPE)[grep("PBMC",names(AMPE))]
head(AMPE)
names(AMPE)

#View(d.anno)

write.csv(d.anno,"csv/d.anno.csv")
head(d.anno)
names(d.anno) <- gsub('.x','',names(d.anno))
head(d.anno)

