library(dplyr)
library(Seurat)
library(ggplot2)

output_folder <- "/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/"

## Must add findvariablefeatures before running as.loom
## Can add entirety of seurat datasets in future

input_folder = "/Users/parisasamareh/Desktop/mount/single_cell_outputs/DAY0_outs/outs/filtered_feature_bc_matrix"
CART_exh.data <- Read10X(data.dir = input_folder)
CART_exh <- CreateSeuratObject(counts = CART_exh.data, project = "CART_exh_Day0", min.cells = 3, min.features = 200)
CART_exh <- FindVariableFeatures(CART_exh)
as.loom(CART_exh, filename = paste0(output_folder, "Day0.loom"))
CART_exh[["percent.mt"]] <- PercentageFeatureSet(CART_exh, pattern = "^MT-")
CART_exh <- subset(CART_exh, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
as.loom(CART_exh, filename = paste0(output_folder, "Day0_filtered.loom"))

rm(list=ls()) 

output_folder <- "/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/"

input_folder = "/Users/parisasamareh/Desktop/mount/single_cell_outputs/DAY20_outs/outs/filtered_feature_bc_matrix"
CART_exh.data <- Read10X(data.dir = input_folder)
CART_exh <- CreateSeuratObject(counts = CART_exh.data, project = "CART_exh_Day20", min.cells = 3, min.features = 200)
CART_exh <- FindVariableFeatures(CART_exh)
as.loom(CART_exh, filename = paste0(output_folder, "Day20.loom"))
CART_exh[["percent.mt"]] <- PercentageFeatureSet(CART_exh, pattern = "^MT-")
CART_exh <- subset(CART_exh, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
as.loom(CART_exh, filename = paste0(output_folder, "Day20_filtered.loom"))

#------------------------------------------------------------------------------------------------------------------------
# Run the jupyter notebook julia file
#------------------------------------------------------------------------------------------------------------------------

rm(list=ls()) 
filtered = TRUE
k = 10
if (filtered && k == 10){
  curr <- "filtered_10"} else if(filtered && k == 20)
    {curr <- "filtered_20"} else if(!filtered && k==10)
    {curr <- "10"} else if(!filtered && k==20)
    {curr <- "20"}  

## Read in files
diff_names = read.table(paste("/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/diff_names_", ".tsv", sep = curr))
diff_pos= read.table(paste("/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/diff_pos_", ".tsv", sep = curr))
diff_neg= read.table(paste("/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/diff_neg_", ".tsv", sep = curr))

## Genes upregulated
diff_pos_ind <- 1*(diff_pos < -1.3)
genes_pos <- rowSums(diff_pos_ind)
## Genes downregulated
diff_neg_ind <- 1*(diff_neg < -1.3)
genes_neg <- rowSums(diff_neg_ind)

## Write df to file
df_genes <- data.frame(genes = diff_names[,1], pos = genes_pos, neg = genes_neg)
write.table(df_genes, 
            file = paste("/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/df_genes_", ".tsv", sep = curr), 
            sep = "\t", quote = FALSE, row.names = FALSE)

## Histogram of upreg
p1 <- ggplot(data=head(df_genes[order(df_genes$pos, decreasing = TRUE),], 100),  aes(x = reorder(genes, -pos), y = pos) )+geom_bar(stat="identity")+ 
  ggtitle("# of Cells Pos. DE vs Gene") + ylab("# of Cells") + xlab("Gene")+ theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))
ggsave(paste("/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/num_cells_vs_gene_pos_", ".pdf", sep = curr), p1)

## Hist of downreg
p2 <- ggplot(data=head(df_genes[order(df_genes$neg, decreasing = TRUE),] , 100),  aes(x = reorder(genes, -neg), y = neg) )+geom_bar(stat="identity")+ 
  ggtitle("# of Cells Neg. DE vs Gene") + ylab("# of Cells") + xlab("Gene")+ theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))
ggsave(paste("/Users/parisasamareh/Desktop/mount/CART_exh_sc_analysis/cellfishing_output/num_cells_vs_gene_neg_", ".pdf", sep = curr))

