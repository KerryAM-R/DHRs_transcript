# Get log2 counts per million
rownames(y.all$genes) <- make.unique(y.all$genes$Gene.Name)
head(y.all$genes)
logcounts <- cpm(y.all,log=F)
rownames(logcounts) <- make.unique(y.all$genes$Gene.Name)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")



# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
names(var_genes) 
var_genes

RSJS <- read.csv("csv/tol.SJS.CD3.trended.csv")

RSJS[abs(RSJS$logFC) > 0.58,]

RSJS.sig <- subset(RSJS,RSJS$Pvalue<0.1 & abs(RSJS$logFC) > 0.58)
dim(RSJS.sig) - 498

RMPE <- read.csv("csv/tol.rMPE.CD3.trended.csv")
RMPE.sig <- subset(RMPE,RMPE$Pvalue<0.1  & abs(RMPE$logFC) > 0.58)
dim(RMPE.sig) -70


AMPE <- read.csv("csv/tol.aMPE.CD3.trended.csv")
AMPE.sig <- subset(AMPE,AMPE$Pvalue<0.1  & abs(AMPE$logFC) > 0.58)
dim(AMPE.sig) - 158

sig.all <- rbind(RSJS.sig,RMPE.sig,AMPE.sig)

unique(sig.all$ID)
library(pheatmap)
library(RColorBrewer)
length(select_var)

var_genes <- apply(logcounts, 1, var)

select_var <- names(sort(var_genes, decreasing=TRUE))[names(sort(var_genes, decreasing=TRUE)) %in% unique(sig.all$ID)]

highly_variable_lcpm <- logcounts[select_var,]

colnames(highly_variable_lcpm) <- paste0(names.all.2$V1,"_",names.all.2$V2)

# Rename columns based on group
colnames(expr_matrix) <- paste0(group_labels, "_", seq(1, 20))
expr_matrix
# Use a color palette for better visualization
colors <- colorRampPalette(brewer.pal(9, "Blues"))(10)
colors
# Ensure there are no missing values in your annotation data
annotation_col <- data.frame(Group = names.all.2$V1)
rownames(annotation_col) <- colnames(highly_variable_lcpm) # Make sure rownames match the column names of the matrix
highly_variable_lcpm

group_colors
annotation_col
dev.off()
# Define colors for each group
group_colors <- c("Tol" = "#BEBEBE", "RMPE" = "#0076C0", "AMPE" = "#E37C1D", "RSJS" = "#A30234")
annotation_colors <- list(Group = group_colors)

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)

# Scale the rows (genes) of the logcounts matrix
scaled_lcpm <- t(scale(t(highly_variable_lcpm)))  # Scale by row

# Define a color function for the heatmap
color_fun <- colorRamp2(c(min(scaled_lcpm), 0, max(scaled_lcpm)), 
                        c("blue", "white", "red"))

# Define the column group annotations
ha_column <- HeatmapAnnotation(Group = names.all.2$V1, 
                               col = list(Group = group_colors),
                               annotation_legend_param = list(
                                 title = "Groups", 
                                 ncol = 1,  # One column for the group legend
                                 legend_direction = "vertical"  # Stack vertically
                               ))

# Create the heatmap WITH row scaling
ht <- Heatmap(scaled_lcpm, 
              name = "Expression", 
              col = color_fun, 
              cluster_rows = TRUE, 
              cluster_columns = T, 
              show_row_names = FALSE, 
              show_column_names = TRUE, 
              column_title = "Samples",
              top_annotation = ha_column,
              heatmap_legend_param = list(
                title = "Scaled Expression", 
                legend_direction = "vertical"  # Vertical direction for heatmap legend
              )
)

# Draw the heatmap and merge legends into one column
draw(ht, 
     merge_legend = TRUE)  # Stack both the annotation and heatmap legends in one column
