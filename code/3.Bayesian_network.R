library(bnlearn)
library(dplyr)
library(igraph)
library(stringr)
library(visNetwork)

Cluster_Dat <- read.csv("../data/curated_pain_data.csv", header = T, stringsAsFactors = T)
df_cluster <- Cluster_Dat %>% mutate_if(is.numeric,as.factor)
all_nodes <- colnames(Cluster_Dat)

# Hill-climbing to learn the network structure
dag.test = hc(df_cluster, restart = 5, perturb = 5, score = 'bic')
# Fit the learned structure to the data
fitted = bn.fit(dag.test, df_cluster, method = "bayes")
# Moralize
moral_net <- moral(dag.test)

# Transform the bn.fit object to igraph object for visualization
g <- igraph.from.graphNEL(as.graphNEL(moral_net))
data <- toVisNetworkData(g)
df_edges <- data$edges

# To moralize edges:
combined <- c()
for(i in 1:length(df_edges$from)){
  combined <- c(combined, paste(sort(c(df_edges[i,'from'], df_edges[i,'to'])), collapse=','))  
}
df_edges$combined <- combined
df_edges <- df_edges %>%
  group_by(combined) %>%
  mutate(id = row_number()) %>%
  ungroup %>%
  subset(id == 1) %>%
  select(to, from, weight)

nod_net <- as.data.frame(data$nodes)
node_labels <- read.csv('display_names.csv')
nod_net <- nod_net %>%
  inner_join(node_labels, by = c('id' = 'column')) %>%
  select(id, display_name) %>%
  rename(label = display_name)

network <- visNetwork(nod_net, edges = df_edges, idToLabel=T,
                     width = "100%", height = 900) %>%
  visIgraphLayout(layout = 'layout_with_kk', type = "full") %>%
  visInteraction(navigationButtons = TRUE,
                 dragNodes = TRUE,
                 dragView = FALSE,
                 zoomView = FALSE,
                 hover=TRUE,
                 hoverConnectedEdges = TRUE) %>%
  visEdges(smooth = list(enabled = TRUE, type = "diagonalCross",roundness = 0.1),
           physics = F,
           color = list(highlight = "blue", hover = "blue")) %>%
  visNodes(color = list(background = "azure",
                        border = "grey",
                        highlight = "firebrick",
                        hover = 'firebrick')) %>%
  visOptions(highlightNearest = list(enabled =TRUE, degree = 1,hover=F,labelOnly=T),
             nodesIdSelection = TRUE)

# visSave(network, file = "../results/bayes_net/Bayes_net_moralized.html", background = "#F5F4F4")
