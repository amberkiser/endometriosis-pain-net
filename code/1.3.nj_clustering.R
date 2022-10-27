# Cluster pain map areas using the neighbor-joining algorithm and Jaccard distance metric.

library(DECIPHER)

# read in CSV with distance matrix (Jaccard distance)
distances <- read.csv('jaccard_distance_painful_areas.csv', row.names=1)
distance_matrix <- data.matrix(distances)

clusters <- IdClusters(distance_matrix, method='NJ', type='clusters', cutoff=0.6375)

clusters <- cbind(rownames(clusters), clusters)
rownames(clusters) <- NULL
colnames(clusters) <- c('painful_area', 'cluster')
write.csv(clusters, '../results/clustering/nj_clustering.csv', row.names = FALSE)
