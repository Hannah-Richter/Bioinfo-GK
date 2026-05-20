library(tidyverse)
library(dplyr)
library(ggplot2)

load("LUAD_data.RData")

expression_data = luad_exp_clean


#qq-plots for all genes from one patient
ggplot(expression_data) +
  geom_qq( aes(sample = get(colnames(expression_data)[2]))) +
  geom_qq( aes(sample = get(colnames(expression_data)[21]))) +
  geom_qq( aes(sample = get(colnames(expression_data)[85]))) +
  geom_qq( aes(sample = get(colnames(expression_data)[107]))) +
  geom_qq( aes(sample = get(colnames(expression_data)[171])))


#Transpose the dataset to allow for gene-wise qq-plots
t_expression_data= t(column_to_rownames(expression_data, "gene")) %>% as.data.frame() %>% rownames_to_column(.,"patient") %>% tibble()
ggplot(t_expression_data) +
  geom_qq(aes(sample = get(colnames(t_expression_data)[4])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[10])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[400])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[566])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[1020])))

#Histograms of some genes
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[6])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[15])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[471])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[1132])))


#Shapiro-Wilks test of all genes separately
norm_test = expression_data[-1] %>% apply(MARGIN = 1, function(x) shapiro.test(x))

#Extract p-values from "Shapiro Test" Object and add to list; plot into histogram
Shapiro_p_values = list()
for(i in 1:length(norm_test)) Shapiro_p_values[i]= norm_test[[i]]$p.value
Shapiro_p_values %>% unlist(use.names = FALSE) %>% hist(,breaks=20)


#Percentage of Shapiro-Wilk Tests whith p smaller 0.05
sum(Shapiro_p_values < 0.05)/length(Shapiro_p_values)

#Quantile Normalisation

#Compute ranks of genes for each patient
ranks <- apply(expression_data[-1], 2, rank, ties.method = "min")

#Sort logCPM values of each patients and calculate average 
sorted_expression = apply(expression_data[-1],2,sort)
#Calculate Row Averages
ranked_means <- rowMeans(sorted_expression)
#Replace Ranks of genes with the values of the means of ranks
norm_expression = as.data.frame(matrix(ranked_means[ranks], ncol=ncol(ranks)))
#Add row- and colnames
rownames(norm_expression) = expression_data[[1]]
colnames(norm_expression) = colnames(expression_data[-1])

#Compare Shapiro-Wilks test
norm_test_normalised = norm_expression %>% apply(MARGIN = 1, function(x) shapiro.test(x))
Shapiro_p_values_normalised = list()
for(i in 1:length(norm_test_normalised)) Shapiro_p_values_normalised[i]= norm_test_normalised[[i]]$p.value
Shapiro_p_values_normalised %>% unlist(use.names = FALSE) %>% hist(,breaks=20)
sum(Shapiro_p_values_normalised < 0.05)/length(Shapiro_p_values_normalised)
#Less normally distributed than before


#Differential Gene Expression

#Which Samples have which subtype?
Prox_inflam = which(luad_anot_clean$paper_expression_subtype=="prox.-inflam")
Prox_prolif = which(luad_anot_clean$paper_expression_subtype=="prox.-prolif.")
Terminal_resp = which(luad_anot_clean$paper_expression_subtype=="TRU")


Exp_inflam = expression_data[Prox_inflam+1]
Exp_prolif = expression_data[Prox_prolif+1]
Exp_term = expression_data[Terminal_resp+1]

# 
sum(!is.na(luad_anot_clean$paper_expression_subtype))
sum(substr(luad_anot_clean$barcode, 14, 15)=="01"&!is.na(luad_anot_clean$paper_expression_subtype))


