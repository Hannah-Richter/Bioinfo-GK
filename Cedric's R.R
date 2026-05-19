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
  geom_qq( aes(sample = get(colnames(expression_data)[171]))) +
  geom_qq( aes(sample = get(colnames(expression_data)[298])))


#Transpose the dataset to allow for gene-wise qq-plots
t_expression_data= t(column_to_rownames(expression_data, "gene")) %>% as.data.frame() %>% rownames_to_column(.,"patient") %>% tibble()
ggplot(t_expression_data) +
  geom_qq(aes(sample = get(colnames(t_expression_data)[4])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[10])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[400])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[566])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[1020])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[2193])))

#Histograms of some genes
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[6])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[15])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[471])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[1132])))
ggplot(t_expression_data) +
  geom_histogram(aes(x = get(colnames(t_expression_data)[2378])))


#Shapiro-Wilks test of all genes separately
norm_test = expression_data[,seq(2,length(expression_data[1,]))] %>% apply(MARGIN =1, function(x) shapiro.test(x))

#Extract p-values from "Shapiro Test" Object and add to list; plot into histogram
Shapiro_p_values = list()
for(i in 1:length(norm_test)) Shapiro_p_values[i]= norm_test[[i]]$p.value
Shapiro_p_values %>% unlist(use.names = FALSE) %>% hist(,breaks=20)


#Percentage of Shapiro-Wilk Tests whith p smaller 0.05
sum(Shapiro_p_values < 0.05)/length(Shapiro_p_values)


rnk <- apply(expression_data[-1], 2, rank, ties.method = "min")
view(rnk)

sorted_expression = apply(expression_data[-1],2,sort)
view(sorted_expression)
ranked_means <- rowMeans(sorted_expression)
view(ranked_means)

norm_expression = as.data.frame(matrix(ranked_means[rnk], ncol=ncol(rnk)))
view(norm_expression)

ggplot(norm_expression) +
  geom_qq( aes(sample = get(colnames(norm_expression)[2])))

t_norm_exp = t(norm_expression) %>% data.frame()

ggplot(t_norm_exp) +
  geom_qq( aes(sample = get(colnames(t_norm_exp)[2])))
