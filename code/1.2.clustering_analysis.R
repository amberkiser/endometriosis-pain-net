# Build chart with number of clusters on X axis and "threshold" used to create those clusters on y-axis.

library(DECIPHER)
library(cluster)
library(dplyr)
library(ggplot2)

distances <- read.csv('jaccard_distance_painful_areas.csv', row.names=1)
distance_matrix <- data.matrix(distances)

cutoffs <- seq(from=0.5, to=1, by=0.0001)
num_of_clusters <- c()
sil_scores <- c()

for(i in 1:length(cutoffs)) {
  clust <- IdClusters(distance_matrix, method='NJ', type='clusters', cutoff=cutoffs[i], verbose=FALSE)
  num_of_clusters <- c(num_of_clusters, max(clust$cluster))
}

cluster_metrics <- data.frame(threshold = cutoffs,
                              number_of_clusters = num_of_clusters)

cluster_number_data <- cluster_metrics %>%
  group_by(number_of_clusters) %>%
  arrange(desc(threshold)) %>%
  mutate(row_number = row_number()) %>%
  filter(row_number == 1) %>%
  select(-row_number) %>%
  ungroup()

ggplot(cluster_number_data, aes(x=number_of_clusters, y=threshold)) +
  geom_line() +
  theme_bw() +
  geom_vline(xintercept = 15, linetype='dashed') +
  labs(x='Number of Clusters', y='Threshold')

ggsave('../results/clustering/cluster_elbow_plot.png')
