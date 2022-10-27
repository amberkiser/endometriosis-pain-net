# Bootstrapping 

library(bnlearn)
library(dplyr)
library(stringr)
library(gRain)
library(xlsx)

Cluster_Dat <- read.csv("../data/curated_pain_data.csv", header = T, stringsAsFactors = T)
Cluster_Dat['index'] <- 1:473
##Convert Numeric to Factor
df_cluster <- Cluster_Dat %>% mutate_if(is.numeric,as.factor)
all_nodes <- colnames(Cluster_Dat)

# Start iteration-----------------------------------------------------------------------------
N <- 1000
probabilities <- data.frame()

for (boot_iter in 1:N) {
  # Bootstrap data
  boot_data <- sample_n(df_cluster, 473, replace=TRUE)
  ib_index <- as.numeric(unique(sort(boot_data$index)))
  oob_index <- unique(Cluster_Dat$index[! Cluster_Dat$index %in% ib_index])
  boot_data <- subset(boot_data, select = -index)
  
  # Hill-climbing to learn the network structure
  dag.test <- hc(boot_data, restart = 5, perturb = 5, score = 'bic')
  
  # Fit the learned structure to the data
  fitted <- bn.fit(dag.test, boot_data, method = "bayes")
  PanelAjunction_AN2 <- compile(as.grain(fitted), propagate = T)

  # Probabilities and relative probabilities
  dx_list <- c('0', '1', '2', '3', '4')
  names(dx_list) <- c('normal', 'endo', 'fibroids', 'cysts', 'other')
  probs <- data.frame()
  
  for (i in 1:length(dx_list)) {
    yes_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c('final_diagnosis_DV'), states = c(dx_list[i])),
                           nodes = "subfertility")$subfertility[2]
    # if dx is normal, do something different
    if(dx_list[i] == '0'){
      no_prob <- 1-yes_prob 
    } else {
      no_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c('final_diagnosis_DV'), states = c('0')),
                            nodes = "subfertility")$subfertility[2]  
    }
    
    # risk of subfertility given dx / risk of subfertility given no dx
    rr <- yes_prob/no_prob

    probs <- rbind(probs, data.frame(diagnosis=names(dx_list[i]), rr=rr))}
  
  probs['iteration'] <- boot_iter
  probabilities <- rbind(probabilities, probs)
} # end of bootstrapped iteration

write.csv(probabilities, '../results/sensitivity_analyses/bootstrap_subfertility.csv', row.names = FALSE)
