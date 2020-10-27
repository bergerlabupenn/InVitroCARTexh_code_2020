library(igraph)
library(pheatmap)

## First create adjacency matrix from network, then visualize using a heatmap
# Read in network file created from PIDC
adj_mat <- read.csv("/Users/parisasamareh/Desktop/temp_gene_exp/day20_network.csv", header=F, as.is=T, sep = "\t")
# Create list of all of the nodes (genes)
nodes <- unique(c(adj_mat$V1, adj_mat$V2))
# Save the network as an igraph object
net <- graph.data.frame(adj_mat[c("V1", "V2")], nodes, directed=F)
# Set the weights of the edges as the strength of the relationship between two genes
E(net)$weight <- adj_mat$V3
# Remove duplicate edges
g <- simplify(net)
# Delete edges with weight = 0 (no relationship between genes)
g <- delete.edges(g, which(E(g)$weight==0))
# Convert the graph to an adjacency matrix
g_adj <- get.adjacency(g, attr="weight", sparse=F)
# Write the adjacency matrix to a csv for future use
write.table(g_adj, paste("/Users/parisasamareh/Desktop/temp_gene_exp", "/day20_adj_mat.csv", sep = ""), sep = "\t")
# Plot the adjacency matrix as a heatmap 
pheatmap(g_adj, fontsize_row = 1.5, fontsize_col = 1.5)