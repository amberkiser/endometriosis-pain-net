library(stats)
library(Exact)
library(xlsx)

# read in data
pain_data <- read.csv('../data/curated_pain_data.csv')

# All Pain ---------------------------------------------------------------------------------------
endo <- pain_data[pain_data$final_diagnosis_DV == 1,]
normal <- pain_data[pain_data$final_diagnosis_DV == 0, ]
fibroids <- pain_data[pain_data$final_diagnosis_DV == 2, ]
cysts <- pain_data[pain_data$final_diagnosis_DV == 3, ]
other <- pain_data[pain_data$final_diagnosis_DV == 4, ]

total_endo <- length(endo$final_diagnosis_DV)
total_normal <- length(normal$final_diagnosis_DV)
total_fibroids <- length(fibroids$final_diagnosis_DV)
total_cysts <- length(cysts$final_diagnosis_DV)
total_other <- length(other$final_diagnosis_DV)
nodes_to_test <- colnames(pain_data)[colnames(pain_data) != 'final_diagnosis_DV']
N <- length(nodes_to_test)


### Fisher Exact test
fisher_multi_df <- data.frame(variable=rep(NA, N), number_endo=rep(NA, N), number_normal=rep(NA, N),
                              number_fibroids=rep(NA, N), number_cysts=rep(NA, N),
                              number_other=rep(NA, N),prop_endo=rep(NA, N), prop_normal=rep(NA, N),
                              prop_fibroids=rep(NA, N), prop_cysts=rep(NA, N),
                              prop_other=rep(NA, N), p_value=rep(NA, N), p_value_adj=rep(NA, N),
                              stringsAsFactors=FALSE) 

for(i in 1:length(nodes_to_test)){  
  endo_with_pain <- length(endo[endo[nodes_to_test[i]] == 1, nodes_to_test[i]])
  endo_no_pain <- length(endo[endo[nodes_to_test[i]] == 0, nodes_to_test[i]])
  normal_with_pain <- length(normal[normal[nodes_to_test[i]] == 1, nodes_to_test[i]])
  normal_no_pain <- length(normal[normal[nodes_to_test[i]] == 0, nodes_to_test[i]])
  fibroids_with_pain <- length(fibroids[fibroids[nodes_to_test[i]] == 1, nodes_to_test[i]])
  fibroids_no_pain <- length(fibroids[fibroids[nodes_to_test[i]] == 0, nodes_to_test[i]])
  cysts_with_pain <- length(cysts[cysts[nodes_to_test[i]] == 1, nodes_to_test[i]])
  cysts_no_pain <- length(cysts[cysts[nodes_to_test[i]] == 0, nodes_to_test[i]])
  other_with_pain <- length(other[other[nodes_to_test[i]] == 1, nodes_to_test[i]])
  other_no_pain <- length(other[other[nodes_to_test[i]] ==0, nodes_to_test[i]])
  contingency_table <- cbind(c(endo_with_pain, endo_no_pain), c(normal_with_pain, normal_no_pain),
                             c(fibroids_with_pain, fibroids_no_pain), c(cysts_with_pain, cysts_no_pain),
                             c(other_with_pain, other_no_pain))
  
  p_val <- fisher.test(contingency_table, workspace=2e7)$p.value
  fisher_multi_df[i,] <- list(nodes_to_test[i], endo_with_pain, normal_with_pain, 
                              fibroids_with_pain, cysts_with_pain, other_with_pain, 
                              endo_with_pain/total_endo, normal_with_pain/total_normal,
                              fibroids_with_pain/total_fibroids, cysts_with_pain/total_cysts,
                              other_with_pain/total_other, p_val, -1)
}
fisher_multi_df$p_value_adj <- p.adjust(fisher_multi_df$p_value, method = 'fdr')

#### Boschloo pairwise comparisons for significant clusters---------------------------------------------------
sig_nodes <- fisher_multi_df[fisher_multi_df$p_value_adj < 0.05, 'variable']
N <- length(sig_nodes)
#### endo vs normal
boschloo_endo_normal_df <- data.frame(variable=rep(NA, N), number_endo=rep(NA, N), 
                                      number_normal=rep(NA, N),prop_endo=rep(NA, N), 
                                      prop_normal=rep(NA, N), p_value=rep(NA, N),  
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE) 
for(i in 1:length(sig_nodes)){
  endo_with_pain <- length(endo[endo[sig_nodes[i]] == 1, sig_nodes[i]])
  endo_no_pain <- length(endo[endo[sig_nodes[i]] == 0, sig_nodes[i]])
  normal_with_pain <- length(normal[normal[sig_nodes[i]] == 1, sig_nodes[i]])
  normal_no_pain <- length(normal[normal[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(endo_with_pain, endo_no_pain), c(normal_with_pain, normal_no_pain))
  
  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo", 
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)
  
  boschloo_endo_normal_df[i,] <- list(sig_nodes[i], endo_with_pain, normal_with_pain, 
                                       endo_with_pain/total_endo, normal_with_pain/total_normal, 
                                       b_test$p.value, -1)
}
boschloo_endo_normal_df$p_value_adj <- p.adjust(boschloo_endo_normal_df$p_value, method = 'fdr')

#### endo vs fibroids
boschloo_endo_fibroids_df <- data.frame(variable=rep(NA, N), number_endo=rep(NA, N),
                                      number_fibroids=rep(NA, N),prop_endo=rep(NA, N),
                                      prop_fibroids=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  endo_with_pain <- length(endo[endo[sig_nodes[i]] == 1, sig_nodes[i]])
  endo_no_pain <- length(endo[endo[sig_nodes[i]] == 0, sig_nodes[i]])
  fibroids_with_pain <- length(fibroids[fibroids[sig_nodes[i]] == 1, sig_nodes[i]])
  fibroids_no_pain <- length(fibroids[fibroids[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(endo_with_pain, endo_no_pain), c(fibroids_with_pain, fibroids_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_endo_fibroids_df[i,] <- list(sig_nodes[i], endo_with_pain, fibroids_with_pain,
                                      endo_with_pain/total_endo, fibroids_with_pain/total_fibroids,
                                      b_test$p.value, -1)
}
boschloo_endo_fibroids_df$p_value_adj <- p.adjust(boschloo_endo_fibroids_df$p_value, method = 'fdr')

#### endo vs cysts
boschloo_endo_cysts_df <- data.frame(variable=rep(NA, N), number_endo=rep(NA, N),
                                      number_cysts=rep(NA, N),prop_endo=rep(NA, N),
                                      prop_cysts=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  endo_with_pain <- length(endo[endo[sig_nodes[i]] == 1, sig_nodes[i]])
  endo_no_pain <- length(endo[endo[sig_nodes[i]] == 0, sig_nodes[i]])
  cysts_with_pain <- length(cysts[cysts[sig_nodes[i]] == 1, sig_nodes[i]])
  cysts_no_pain <- length(cysts[cysts[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(endo_with_pain, endo_no_pain), c(cysts_with_pain, cysts_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_endo_cysts_df[i,] <- list(sig_nodes[i], endo_with_pain, cysts_with_pain,
                                      endo_with_pain/total_endo, cysts_with_pain/total_cysts,
                                      b_test$p.value, -1)
}
boschloo_endo_cysts_df$p_value_adj <- p.adjust(boschloo_endo_cysts_df$p_value, method = 'fdr')

#### endo vs other
boschloo_endo_other_df <- data.frame(variable=rep(NA, N), number_endo=rep(NA, N), 
                                     number_other=rep(NA, N),prop_endo=rep(NA, N), 
                                     prop_other=rep(NA, N), p_value=rep(NA, N),
                                     p_value_adj=rep(NA, N), stringsAsFactors=FALSE) 
for(i in 1:length(sig_nodes)){
  endo_with_pain <- length(endo[endo[sig_nodes[i]] == 1, sig_nodes[i]])
  endo_no_pain <- length(endo[endo[sig_nodes[i]] == 0, sig_nodes[i]])
  other_with_pain <- length(other[other[sig_nodes[i]] == 1, sig_nodes[i]])
  other_no_pain <- length(other[other[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(endo_with_pain, endo_no_pain), c(other_with_pain, other_no_pain))
  
  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo", 
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)
  
  boschloo_endo_other_df[i,] <- list(sig_nodes[i], endo_with_pain, other_with_pain, 
                                     endo_with_pain/total_endo, other_with_pain/total_other, 
                                     b_test$p.value, -1)
}
boschloo_endo_other_df$p_value_adj <- p.adjust(boschloo_endo_other_df$p_value, method = 'fdr')

#### fibroids vs cysts
boschloo_fibroids_cysts_df <- data.frame(variable=rep(NA, N), number_fibroids=rep(NA, N),
                                      number_cysts=rep(NA, N),prop_fibroids=rep(NA, N),
                                      prop_cysts=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  fibroids_with_pain <- length(fibroids[fibroids[sig_nodes[i]] == 1, sig_nodes[i]])
  fibroids_no_pain <- length(fibroids[fibroids[sig_nodes[i]] == 0, sig_nodes[i]])
  cysts_with_pain <- length(cysts[cysts[sig_nodes[i]] == 1, sig_nodes[i]])
  cysts_no_pain <- length(cysts[cysts[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(fibroids_with_pain, fibroids_no_pain), c(cysts_with_pain, cysts_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_fibroids_cysts_df[i,] <- list(sig_nodes[i], fibroids_with_pain, cysts_with_pain,
                                      fibroids_with_pain/total_fibroids, cysts_with_pain/total_cysts,
                                      b_test$p.value, -1)
}
boschloo_fibroids_cysts_df$p_value_adj <- p.adjust(boschloo_fibroids_cysts_df$p_value, method = 'fdr')

#### fibroids vs other
boschloo_fibroids_other_df <- data.frame(variable=rep(NA, N), number_fibroids=rep(NA, N),
                                      number_other=rep(NA, N),prop_fibroids=rep(NA, N),
                                      prop_other=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  fibroids_with_pain <- length(fibroids[fibroids[sig_nodes[i]] == 1, sig_nodes[i]])
  fibroids_no_pain <- length(fibroids[fibroids[sig_nodes[i]] == 0, sig_nodes[i]])
  other_with_pain <- length(other[other[sig_nodes[i]] == 1, sig_nodes[i]])
  other_no_pain <- length(other[other[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(fibroids_with_pain, fibroids_no_pain), c(other_with_pain, other_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_fibroids_other_df[i,] <- list(sig_nodes[i], fibroids_with_pain, other_with_pain,
                                      fibroids_with_pain/total_fibroids, other_with_pain/total_other,
                                      b_test$p.value, -1)
}
boschloo_fibroids_other_df$p_value_adj <- p.adjust(boschloo_fibroids_other_df$p_value, method = 'fdr')

#### fibroids vs normal
boschloo_fibroids_normal_df <- data.frame(variable=rep(NA, N), number_fibroids=rep(NA, N),
                                      number_normal=rep(NA, N),prop_fibroids=rep(NA, N),
                                      prop_normal=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  fibroids_with_pain <- length(fibroids[fibroids[sig_nodes[i]] == 1, sig_nodes[i]])
  fibroids_no_pain <- length(fibroids[fibroids[sig_nodes[i]] == 0, sig_nodes[i]])
  normal_with_pain <- length(normal[normal[sig_nodes[i]] == 1, sig_nodes[i]])
  normal_no_pain <- length(normal[normal[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(fibroids_with_pain, fibroids_no_pain), c(normal_with_pain, normal_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_fibroids_normal_df[i,] <- list(sig_nodes[i], fibroids_with_pain, normal_with_pain,
                                      fibroids_with_pain/total_fibroids, normal_with_pain/total_normal,
                                      b_test$p.value, -1)
}
boschloo_fibroids_normal_df$p_value_adj <- p.adjust(boschloo_fibroids_normal_df$p_value, method = 'fdr')

#### cysts vs other
boschloo_cysts_other_df <- data.frame(variable=rep(NA, N), number_cysts=rep(NA, N),
                                      number_other=rep(NA, N),prop_cysts=rep(NA, N),
                                      prop_other=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  cysts_with_pain <- length(cysts[cysts[sig_nodes[i]] == 1, sig_nodes[i]])
  cysts_no_pain <- length(cysts[cysts[sig_nodes[i]] == 0, sig_nodes[i]])
  other_with_pain <- length(other[other[sig_nodes[i]] == 1, sig_nodes[i]])
  other_no_pain <- length(other[other[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(cysts_with_pain, cysts_no_pain), c(other_with_pain, other_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_cysts_other_df[i,] <- list(sig_nodes[i], cysts_with_pain, other_with_pain,
                                      cysts_with_pain/total_cysts, other_with_pain/total_other,
                                      b_test$p.value, -1)
}
boschloo_cysts_other_df$p_value_adj <- p.adjust(boschloo_cysts_other_df$p_value, method = 'fdr')

#### cysts vs normal
boschloo_cysts_normal_df <- data.frame(variable=rep(NA, N), number_cysts=rep(NA, N),
                                      number_normal=rep(NA, N),prop_cysts=rep(NA, N),
                                      prop_normal=rep(NA, N), p_value=rep(NA, N),
                                      p_value_adj=rep(NA, N), stringsAsFactors=FALSE)
for(i in 1:length(sig_nodes)){
  cysts_with_pain <- length(cysts[cysts[sig_nodes[i]] == 1, sig_nodes[i]])
  cysts_no_pain <- length(cysts[cysts[sig_nodes[i]] == 0, sig_nodes[i]])
  normal_with_pain <- length(normal[normal[sig_nodes[i]] == 1, sig_nodes[i]])
  normal_no_pain <- length(normal[normal[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(cysts_with_pain, cysts_no_pain), c(normal_with_pain, normal_no_pain))

  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo",
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)

  boschloo_cysts_normal_df[i,] <- list(sig_nodes[i], cysts_with_pain, normal_with_pain,
                                      cysts_with_pain/total_cysts, normal_with_pain/total_normal,
                                      b_test$p.value, -1)
}
boschloo_cysts_normal_df$p_value_adj <- p.adjust(boschloo_cysts_normal_df$p_value, method = 'fdr')

#### normal vs other
boschloo_normal_other_df <- data.frame(cluster=rep(NA, N), number_normal=rep(NA, N),
                                        number_other=rep(NA, N), prop_normal=rep(NA, N),
                                        prop_other=rep(NA, N), p_value=rep(NA, N),
                                        p_value_adj=rep(NA, N), stringsAsFactors=FALSE) 
for(i in 1:length(sig_nodes)){
  normal_with_pain <- length(normal[normal[sig_nodes[i]] == 1, sig_nodes[i]])
  normal_no_pain <- length(normal[normal[sig_nodes[i]] == 0, sig_nodes[i]])
  other_with_pain <- length(other[other[sig_nodes[i]] == 1, sig_nodes[i]])
  other_no_pain <- length(other[other[sig_nodes[i]] == 0, sig_nodes[i]])
  contingency_table <- cbind(c(normal_with_pain, normal_no_pain), c(other_with_pain, other_no_pain))
  
  b_test <- exact.test(contingency_table, alternative = "two.sided", method = "boschloo", 
                       model = "Binomial", tsmethod = "central", to.plot = FALSE)
  
  boschloo_normal_other_df[i,] <- list(sig_nodes[i], normal_with_pain, other_with_pain, 
                                        normal_with_pain/total_normal, other_with_pain/total_other, 
                                        b_test$p.value, -1)
}
boschloo_normal_other_df$p_value_adj <- p.adjust(boschloo_normal_other_df$p_value, method = 'fdr')

## Save results
## workbook 2 - endo vs healthy vs other
wb2 <- createWorkbook(type="xlsx")

sheet1 <- createSheet(wb2, sheetName = "fisher_exact")
addDataFrame(fisher_multi_df, sheet1, row.names = FALSE)

sheet2 <- createSheet(wb2, sheetName = "boschloo_endoVnormal")
addDataFrame(na.omit(boschloo_endo_normal_df), sheet2, row.names = FALSE)

sheet3 <- createSheet(wb2, sheetName = "boschloo_endoVfibroids")
addDataFrame(na.omit(boschloo_endo_fibroids_df), sheet3, row.names = FALSE)

sheet4 <- createSheet(wb2, sheetName = "boschloo_endoVcysts")
addDataFrame(na.omit(boschloo_endo_cysts_df), sheet4, row.names = FALSE)

sheet5 <- createSheet(wb2, sheetName = "boschloo_endoVother")
addDataFrame(na.omit(boschloo_endo_other_df), sheet5, row.names = FALSE)

sheet6 <- createSheet(wb2, sheetName = "boschloo_fibroidsVcysts")
addDataFrame(na.omit(boschloo_fibroids_cysts_df), sheet6, row.names = FALSE)

sheet7 <- createSheet(wb2, sheetName = "boschloo_fibroidsVother")
addDataFrame(na.omit(boschloo_fibroids_other_df), sheet7, row.names = FALSE)

sheet8 <- createSheet(wb2, sheetName = "boschloo_fibroidsVnormal")
addDataFrame(na.omit(boschloo_fibroids_normal_df), sheet8, row.names = FALSE)

sheet9 <- createSheet(wb2, sheetName = "boschloo_cystsVother")
addDataFrame(na.omit(boschloo_cysts_other_df), sheet9, row.names = FALSE)

sheet10 <- createSheet(wb2, sheetName = "boschloo_cystsVnormal")
addDataFrame(na.omit(boschloo_cysts_normal_df), sheet10, row.names = FALSE)

sheet11 <- createSheet(wb2, sheetName = "boschloo_normalVother")
addDataFrame(na.omit(boschloo_normal_other_df), sheet11, row.names = FALSE)

saveWorkbook(wb2, "../results/sensitivity_analyses/stat_tests.xlsx")
