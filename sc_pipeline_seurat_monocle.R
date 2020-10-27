library(dplyr)
library(Seurat)
library(ggplot2)
library(monocle3)

## Load the datas
# sshfs psamareh@mercury.pmacs.upenn.edu:/home/psamareh ~/Desktop/mount

# Location of counts data -- should have barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
input_folder = "/Users/parisasamareh/Downloads/filtered_gene_bc_matrices"
# Location of the folder you want your plots to be in
output_folder = "~/Desktop/work_for_charly/07_20_2020_metascapeforday0ND388"
# Read in the 10x data
CART_exh.data <- Read10X(data.dir = input_folder)

## Initialize the Seurat object with the raw (non-normalized data).
# Only keep genes expressed in at least 3 cells and cells that express at least 200 genes.
CART_exh <- CreateSeuratObject(counts = CART_exh.data, project = "CART_exh", min.cells = 3, min.features = 200)
CART_exh[["percent.mt"]] <- PercentageFeatureSet(CART_exh, pattern = "^MT-")

## Pre-processing and QC
# Violin plot of number of features per cell, number of counts per cell, 
# and percentage of genes in the cell that map to the mitochondiral genome
pdf(paste(output_folder, "/violin_plot_QC.pdf", sep = ""))
VlnPlot(CART_exh, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
dev.off()

# Correlation between percent mito and RNA count as well as number of genes expressed in a cell and RNA count.
pdf(paste(output_folder, "/feature_scatter_QC.pdf", sep = ""))
plot1 <- FeatureScatter(CART_exh, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CART_exh, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Printing the number of cells we have without the mitochondrial and gene filter
print("Number of cells without mito and gene filter:\n")
print(dim(CART_exh@meta.data))
# Each cell expresses between 200 and 5000 genes and has less than 5% mito dna
CART_exh <- subset(CART_exh, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
# Printing the number of cells we have with the mitochondrial and gene filter
print("Number of cells with mito and gene filter:\n")
print(dim(CART_exh@meta.data))
# Printing out the normalized gene expression data to a csv
write.csv(t(as.matrix(CART_exh[["RNA"]]@data)), paste("/Users/parisasamareh/Desktop/temp_gene_exp", "/gene_exp_day20_mito5.csv", sep = ""))

# Perform SCTransform normalization
CART_exh <- SCTransform(CART_exh)

# Run PCA
CART_exh <- RunPCA(CART_exh, features = VariableFeatures(object = CART_exh))
# If you want to print genes associated with each PC: print(CART_exh[["pca"]], dims = 1:5, nfeatures = 10)

# Plot PC 1 and 2 against each other
pdf(paste(output_folder, "/pca_plot.pdf", sep = ""))
DimPlot(CART_exh, reduction = "pca")
dev.off()

# Plot heatmatmap of up and down regulated genes for each PC using 500 cells
pdf(paste(output_folder, "/heatmap_PCs.pdf", sep = ""))
DimHeatmap(CART_exh, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

# Plot standard deviation each PC encompasses to determine the most informative PCs
pdf(paste(output_folder, "/elbow.pdf", sep = ""))
ElbowPlot(CART_exh)
dev.off()

# Construct a KNN graph based on the euclidean distance in PCA space, 
# and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).
CART_exh <- FindNeighbors(CART_exh, dims = 1:10)
# Apply the Louvain algorithm to cluster cells by iteratively grouping them together, optimizing modularity
# The higher the resolution, the more clusters: "0.4-1.2 typically returns good results for single-cell datasets of around 3K cells"
CART_exh <- FindClusters(CART_exh, resolution = 0.2)
# Run non-linear dimensionality reduction to visualize data: use the same number of PCs as in clustering
CART_exh <- RunUMAP(CART_exh, dims = 1:10)
CART_exh <- RunTSNE(CART_exh, dims = 1:10)

# Plot UMAP
pdf(paste(output_folder, "/umap.pdf", sep = ""))
DimPlot(CART_exh, reduction = "umap")
dev.off()
# Plot tSNE
pdf(paste(output_folder, "/tsne.pdf", sep = ""))
DimPlot(CART_exh, reduction = "tsne")
dev.off()

## Gene markers
# Find the positive markers for each cluster in the single cell data: must be detected in at least 25% of the cells and have a log fc of 0.25
CART_exh.markers <- FindAllMarkers(CART_exh, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Save the markers in a csv 
write.csv(CART_exh.markers, paste(output_folder, "/cluster_markers.csv", sep = ""))
# If you want to find the markers for more than one cluster.
# cluster23.markers <- FindMarkers(CART_exh, ident.1 = c(2,3), min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(cluster23.markers, paste(output_folder, "/cluster23_markers.csv", sep = ""))
CART_exh.markers.03 <- FindMarkers(CART_exh, ident.1 = c(0,3), logfc.threshold = 0)
EnhancedVolcano(CART_exh.markers.03, subtitle = "Differentially Expressed Genes in Clusters 0 and 3",
                lab = rownames(CART_exh.markers.03),
                x = 'avg_logFC',
                y = 'p_val',
                xlim = c(-3, 3), pCutoff = 0.05, FCcutoff = 0.5)

# A shortened list of the top 10 gene markers for each cluster based of log FC
top10 <- CART_exh.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# Heatmap of the gene expression for the top 10 gene markers in each cluster
pdf(paste(output_folder, "/cluster_heatmap.pdf", sep = ""))
DoHeatmap(CART_exh, features = top10$gene) + NoLegend()
dev.off()

## Plots to investigate particular genes
# Create violin plot of expression of these genes, separated by cluster.
pdf(paste(output_folder, "/vlnplot_counts.pdf", sep = ""))
VlnPlot(CART_exh, features = c("GNLY", "ENTPD1","LAYN","KLRC1", "KLRC2", "KLRB1", "KLRD1", "PHLDA1", "SRGAP3"), slot = "counts", log = TRUE)
dev.off()

# UMAP of gene expression for the following genes: one plot saved for each gene
for (i in c("KLRC1", "TOX", 
            "GNLY",  
            "LAYN",
            "CCL4",
            "GNLY",
            "ENTPD1",
            "TIGIT",
            "PHLDA1",
            "GZMA")){
  if(!(i %in% all.genes)){
    print(i)
  }
  pdf(paste(c(output_folder, "/", i, ".pdf"), collapse = ""))
  curr <- FeaturePlot(CART_exh, features = c(i))
  print(curr)
  dev.off()
}

# Dot plot where dot size corresponds to the percentage of cells in that cluster expressing the gene
# and color corresponds to the average gene expression level of those expressing cells
pdf(paste(output_folder, "/dot_plot.pdf", sep = ""))
DotPlot(CART_exh, 
        features = c("TNFRSF18", "GNLY", "ENTPD1", "PHLDA1", 
                      "SRGAP3", "SOX4")) + coord_flip()
dev.off()

#--------------------------------------------------------------------------------------------------------
## Monocle

# Get the raw, subsetted (for mito etc.) counts from the seurat object
data <- as(as.matrix(CART_exh@assays$RNA@counts), 'sparseMatrix')
# Get the cell meta data, mito percentage, seurat cluster, etc.)
cell_metadata <- CART_exh@meta.data
# Gene names
gene_metadata <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
# Create the monocle object
cds <- new_cell_data_set(expression_data = data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# Normalizes data by log and size factor
# then calculates a lower dimensional space using PCA that will be used as input in the future
cds <- preprocess_cds(cds, num_dim = 100)

# Perform dimensionality reduction, UMAP is default
cds <- reduce_dimension(cds)
# If you want to perform tsne as well, cds <- reduce_dimension(cds, reduction_method = 'tSNE')

# Perform clustering: default dim reduction it clusters on is UMAP
cds = cluster_cells(cds)

# Plot UMAP of cells with monocle clusters
pdf(paste(output_folder, "/monocle_umap_clusters.pdf", sep = ""))
plot_cells(cds, reduction_method="UMAP")
dev.off()

# Plot tSNE of cells with monocle clusters
pdf(paste(output_folder, "/monocle_tsne_clusters.pdf", sep = ""))
plot_cells(cds, reduction_method="tSNE")
dev.off()

## Can see how these plots compare to seurat: monocle dim reduction with seurat clusters colored
plot_cells(cds, color_cells_by="seurat_clusters")
plot_cells(cds, color_cells_by="seurat_clusters", reduction_method = "tSNE")
## Plot gene expression on UMAP to see if expression is cluster associated
plot_cells(cds, genes=c("KLRC1"))
plot_cells(cds, genes=c("LAYN"))

## Marker genes for monocle clusters
marker_test_res_mon <- top_markers(cds, group_cells_by="cluster")
write.csv(marker_test_res_mon, file = paste(output_folder, "/monocle_markers_for_monocle_clusters.csv", sep = ""), quote = FALSE)

## Marker genes for seurat clusters
marker_test_res_seur <- top_markers(cds, group_cells_by="seurat_clusters")
write.csv(marker_test_res_seur, file = paste(output_folder, "/monocle_markers_for_seurat_clusters.csv", sep = ""), quote = FALSE)

## Modules for entire data set
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
# Filter genes based on Morans values:
# tells you whether cells at nearby positions on a trajectory will have similar (or dissimilar) expression levels for the gene being tested
# +1 means nearby cells will have perfectly similar expression, -1 is anti-correlated, 0 is no correlation
# Also want siginificant q-value
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))
# Find modules using the genes that showed to be significant with the Morans I test
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution = 0.01)
# Save module csv
write.csv(gene_module_df, file = paste(output_folder, "/modules_all.csv", sep = ""), quote = FALSE)
# Plot umap with module expression for each module
pdf(paste(output_folder, "/modules_umap.pdf", sep = ""))
plot_cells(cds, genes=gene_module_df, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)
dev.off()

## Plot module expression with the seurat UMAP
for (i in levels(gene_module_df$module)){
  seurat_object <- AddModuleScore(seurat_object, list(subset(gene_module_df, module == i)$id), name = paste("mod_", i, sep = ""))
}
x <- colnames(seurat_object@meta.data)

for (i in subset(x, grepl("mod", x))){
  pdf(paste(c(output_folder, "/", i, ".pdf"), collapse = ""))
  curr <- FeaturePlot(seurat_object, features = c(i))
  print(curr)
  dev.off()
}

## Trajectory analysis
cds <- learn_graph(cds)
cds <- order_cells(cds)

pdf(paste(output_folder, "/trajectories_seurat_clusters.pdf", sep = ""))
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE)+ theme(
             legend.position = c(0.95, 0.95),
             legend.justification = c("right", "top")
           )
dev.off()

pdf(paste(output_folder, "/trajectories_monocle_clusters.pdf", sep = ""))
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE)+ theme(
             legend.position = c(0.95, 0.95),
             legend.justification = c("right", "top")
           )
dev.off()

