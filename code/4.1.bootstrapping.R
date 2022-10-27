# Bootstrapping to get risks and relative risks

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
performance_metrics <- data.frame()
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
  
  # Performance metrics
  bic <- score(dag.test, boot_data, type='bic')
  aic <- score(dag.test, boot_data, type='aic')

  X_val <- subset(Cluster_Dat, index %in% oob_index, select = -c(index, final_diagnosis_DV))
  y_val <- subset(Cluster_Dat, index %in% oob_index, select = final_diagnosis_DV)$final_diagnosis_DV
  nodes_to_test <- colnames(X_val)

  y_prob_0 <- c()
  y_prob_1 <- c()
  y_prob_2 <- c()
  y_prob_3 <- c()
  y_prob_4 <- c()
  y_pred <- c()
  for (i in 1:length(X_val$abdominal_pain_12wks)){
    test_states <- as.character(X_val[i,])

    prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = nodes_to_test,
                                  states = test_states),
                       nodes = "final_diagnosis_DV")$final_diagnosis_DV

    y_pred <- c(y_pred, names(which.max(prob)))
    y_prob_0 <- c(y_prob_0, prob[1])
    y_prob_1 <- c(y_prob_1, prob[2])
    y_prob_2 <- c(y_prob_2, prob[3])
    y_prob_3 <- c(y_prob_3, prob[4])
    y_prob_4 <- c(y_prob_4, prob[5])
  }
  val_df <- data.frame(iteration=boot_iter, bic=bic, aic=aic, y_true=y_val, y_pred=y_pred,
                       y_prob_0=y_prob_0, y_prob_1=y_prob_1, y_prob_2=y_prob_2, y_prob_3=y_prob_3, 
                       y_prob_4=y_prob_4)
  performance_metrics <- rbind(performance_metrics, val_df)

  # Probabilities and relative probabilities
  nodes_to_test <- colnames(boot_data)[colnames(boot_data) != 'final_diagnosis_DV']
  # endo---------------------------------------------------------------------------------------
  base_endo <- querygrain(PanelAjunction_AN2, nodes = "final_diagnosis_DV")$final_diagnosis_DV[2]

  endo_probs <- data.frame()
  for (i in 1:length(nodes_to_test)) {
    condition <- nodes_to_test[i]
    yes_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('1')),
                           nodes = "final_diagnosis_DV")$final_diagnosis_DV[2]
    no_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('0')),
                          nodes = "final_diagnosis_DV")$final_diagnosis_DV[2]

    yes_odds <- yes_prob/(1-yes_prob)
    no_odds <- no_prob/(1-no_prob)
    odds_ratio <- yes_odds/no_odds

    rr <- yes_prob/no_prob

    endo_probs <- rbind(endo_probs, data.frame(node=condition, or_endo=odds_ratio, rr_endo=rr))}

  # fibroids---------------------------------------------------------------------------------------
  base_fibroids <- querygrain(PanelAjunction_AN2, nodes = "final_diagnosis_DV")$final_diagnosis_DV[3]
  
  fibroids_probs <- data.frame()
  for (i in 1:length(nodes_to_test)) {
    condition <- nodes_to_test[i]
    yes_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('1')),
                           nodes = "final_diagnosis_DV")$final_diagnosis_DV[3]
    no_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('0')),
                          nodes = "final_diagnosis_DV")$final_diagnosis_DV[3]
    
    yes_odds <- yes_prob/(1-yes_prob)
    no_odds <- no_prob/(1-no_prob)
    odds_ratio <- yes_odds/no_odds
    
    rr <- yes_prob/no_prob
    
    fibroids_probs <- rbind(fibroids_probs, data.frame(node=condition, or_fibroids=odds_ratio, rr_fibroids=rr))}
  
  # cysts---------------------------------------------------------------------------------------
  base_cysts <- querygrain(PanelAjunction_AN2, nodes = "final_diagnosis_DV")$final_diagnosis_DV[4]
  
  cysts_probs <- data.frame()
  for (i in 1:length(nodes_to_test)) {
    condition <- nodes_to_test[i]
    yes_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('1')),
                           nodes = "final_diagnosis_DV")$final_diagnosis_DV[4]
    no_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('0')),
                          nodes = "final_diagnosis_DV")$final_diagnosis_DV[4]
    
    yes_odds <- yes_prob/(1-yes_prob)
    no_odds <- no_prob/(1-no_prob)
    odds_ratio <- yes_odds/no_odds
    
    rr <- yes_prob/no_prob
    
    cysts_probs <- rbind(cysts_probs, data.frame(node=condition, or_cysts=odds_ratio, rr_cysts=rr))}
  
  # other---------------------------------------------------------------------------------------
  base_other <- querygrain(PanelAjunction_AN2, nodes = "final_diagnosis_DV")$final_diagnosis_DV[5]

  other_probs <- data.frame()
  for (i in 1:length(nodes_to_test)) {
    condition <- nodes_to_test[i]
    yes_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('1')),
                           nodes = "final_diagnosis_DV")$final_diagnosis_DV[5]
    no_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('0')),
                          nodes = "final_diagnosis_DV")$final_diagnosis_DV[5]

    yes_odds <- yes_prob/(1-yes_prob)
    no_odds <- no_prob/(1-no_prob)
    odds_ratio <- yes_odds/no_odds

    rr <- yes_prob/no_prob

    other_probs <- rbind(other_probs, data.frame(node=condition, or_other=odds_ratio, rr_other=rr))}

  # normal---------------------------------------------------------------------------------------
  base_normal <- querygrain(PanelAjunction_AN2, nodes = "final_diagnosis_DV")$final_diagnosis_DV[1]

  normal_probs <- data.frame()
  for (i in 1:length(nodes_to_test)) {
    condition <- nodes_to_test[i]
    yes_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('1')),
                           nodes = "final_diagnosis_DV")$final_diagnosis_DV[1]
    no_prob <- querygrain(setFinding(PanelAjunction_AN2, nodes = c(condition), states = c('0')),
                          nodes = "final_diagnosis_DV")$final_diagnosis_DV[1]

    yes_odds <- yes_prob/(1-yes_prob)
    no_odds <- no_prob/(1-no_prob)
    odds_ratio <- yes_odds/no_odds

    rr <- yes_prob/no_prob

    normal_probs <- rbind(normal_probs, data.frame(node=condition, or_normal=odds_ratio, rr_normal=rr))}

  probs <- merge(endo_probs, other_probs, all=TRUE)
  probs <- merge(probs, normal_probs, all=TRUE)
  probs <- merge(probs, fibroids_probs, all=TRUE)
  probs <- merge(probs, cysts_probs, all=TRUE)
  probs['baseline_endo'] <- base_endo
  probs['baseline_other'] <- base_other
  probs['baseline_normal'] <- base_normal
  probs['baseline_fibroids'] <- base_fibroids
  probs['baseline_cysts'] <- base_cysts
  probs['iteration'] <- boot_iter
  probabilities <- rbind(probabilities, probs)
} # end of bootstrapped iteration

write.csv(probabilities, '../results/bayes_net/bootstrap_probabilities.csv', row.names = FALSE)
write.csv(performance_metrics, '../results/bayes_net/bootstrap_predictions.csv', row.names = FALSE)
