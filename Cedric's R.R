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

#Combine dataframes for easier (or maybe actually worse handling) and create list of unique subytpes as well as the subytpe each sample posseses
exp_df = as.data.frame(luad_exp_clean) %>% column_to_rownames("gene") %>% t()
gene_names = colnames(exp_df)
anot_df = as.data.frame(luad_anot_clean) %>% column_to_rownames("barcode")
combined_df = cbind(exp_df,anot_df)
subtypes = unique(anot_df$paper_expression_subtype) %>% na.omit()
sample_subtype = combined_df$paper_expression_subtype


#Prepare results data frame with all the genes and subytpes to be tested vs. rest
DGE_df = data.frame(matrix(NA, nrow = length(gene_names), ncol = length(subtypes)))
rownames(DGE_df) = gene_names
colnames(DGE_df) = subtypes

#Run a Wilcox test for each subtype vs rest across each gene and save p-value of the test into result dataframe
for (st in subtypes){
  not_st = subtypes[subtypes != st]
  for (gene in gene_names) {
  exp_vals = combined_df[[gene]]
  in_group = exp_vals[which(sample_subtype == st)]
  out_group = exp_vals[which(sample_subtype %in% not_st)]
  temporary_test_result = wilcox.test(in_group,out_group,alternative = "two.sided")
  DGE_df[gene,st] = temporary_test_result$p.value
   }
}




