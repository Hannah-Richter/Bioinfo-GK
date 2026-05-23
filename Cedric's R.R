library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)

rm(list=ls())

load("LUAD_data.RData")

expression_data = luad_exp_clean

NA_flattener <- function(data) {
  for (i in seq_along(data)) {
    data[tolower(data[[i]]) %in% c("[unknown]","[not available]"), i ] <- NA
  }
  return(data)
}

luad_anot_clean = NA_flattener(luad_anot_clean)

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
Shapiro_p_values = vector()
for(i in 1:length(norm_test)) Shapiro_p_values[i]= norm_test[[i]]$p.value
Shapiro_p_values %>% hist(,breaks=20)


low_norm_genes = which(rank(Shapiro_p_values)<=5)
medium_norm_genes = which(rank(Shapiro_p_values) %in% seq(1500,1504))
high_norm_genes = which(rank(Shapiro_p_values)>{max(rank(Shapiro_p_values))-5})


t_expression_data= t(column_to_rownames(expression_data, "gene")) %>% as.data.frame() %>% rownames_to_column(.,"patient") %>% tibble()

ggplot(t_expression_data) +
  geom_qq(aes(sample = get(colnames(t_expression_data)[high_norm_genes[1]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[high_norm_genes[2]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[high_norm_genes[3]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[high_norm_genes[4]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[high_norm_genes[5]])))

ggplot(t_expression_data) +
  geom_qq(aes(sample = get(colnames(t_expression_data)[medium_norm_genes[1]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[medium_norm_genes[2]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[medium_norm_genes[3]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[medium_norm_genes[4]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[medium_norm_genes[5]])))

ggplot(t_expression_data) +
  geom_qq(aes(sample = get(colnames(t_expression_data)[low_norm_genes[1]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[low_norm_genes[2]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[low_norm_genes[3]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[low_norm_genes[4]])))+
  geom_qq(aes(sample = get(colnames(t_expression_data)[low_norm_genes[5]])))

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

ggplot({t(norm_expression)}) +
  geom_qq(aes(sample = get(colnames({t(norm_expression)})[low_norm_genes[1]])))+
  geom_qq(aes(sample = get(colnames({t(norm_expression)})[low_norm_genes[2]])))+
  geom_qq(aes(sample = get(colnames({t(norm_expression)})[low_norm_genes[3]])))+
  geom_qq(aes(sample = get(colnames({t(norm_expression)})[low_norm_genes[4]])))+
  geom_qq(aes(sample = get(colnames({t(norm_expression)})[low_norm_genes[5]])))


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

#BH correciton
BH_DGE_df = flatten(DGE_df) %>% p.adjust("BH") %>% matrix(ncol = 3) %>% as.data.frame(row.names = gene_names)
colnames(BH_DGE_df) = subtypes


#Top 18 genes for each subtype vs rest comparison
top_DEGs = gene_names[which(apply(BH_DGE_df,2,rank)<18)%%length(gene_names)]
n_distinct(top_DEGs)
# there are only 49 unique entries among the 51 total entries
#The two repeated genes are:
unique(top_DEGs) %>% sort() %>% .[which(as.vector(table(top_DEGs))>1)]
#Both of these are repeated twice. We will simply make a heatmap of the 49 unique DEGs among the 3 top 18 lists.

#Select top 49 genes
DEG_st_exp = as.data.frame(expression_data) %>% column_to_rownames(colnames(expression_data)[1]) %>% .[expression_data[[1]] %in% unique(top_DEGs),which(!is.na(luad_anot_clean$paper_expression_subtype))]
#Select desired annotation columns
DEG_anot = as.data.frame(luad_anot_clean) %>% column_to_rownames(colnames(luad_anot_clean)[1]) %>% select(paper_expression_subtype,ajcc_pathologic_stage,paper_Smoking.Status) %>% .[which(!is.na(luad_anot_clean$paper_expression_subtype)),]
#Heatmap with annotation
pheatmap(t(DEG_st_exp), show_colnames = FALSE, show_rownames = FALSE, annotation_row = DEG_anot)

#Select top 5 DEGs across all subtypes
BH_DGE_df = BH_DGE_df %>% mutate(minimum = apply(BH_DGE_df, MARGIN = 1, FUN = function(x) min(x)))
#Select only those 5 genes and only the samples for which we know expression subtypes
top_5_subtype = exp_df[,which(rank(BH_DGE_df$minimum)<6)] %>% as.data.frame() %>% mutate(subtype=luad_anot_clean$paper_expression_subtype)
#top_5_expression = exp_df[which(!is.na(luad_anot_clean$paper_expression_subtype)),which(rank(BH_DGE_df$minimum)<6)] %>% as.data.frame() %>% mutate(subtype=DEG_anot$paper_expression_subtype)
top_5_expression_long = top_5_subtype %>% rownames_to_column(var = "Sample ID") %>% pivot_longer(cols = c(colnames(top_5_subtype)[1:5]), names_to = "gene", values_to = "expr")

ggplot(top_5_expression_long, aes(x = subtype, y = expr, fill = subtype)) +
  geom_boxplot() +
  facet_wrap(~gene) +
  labs(x = "Subtype", y = "Expression") +
  theme_bw()

#Oncogenes or TUmour Suppressor
#ADH1B: Alcohol dehydrogenase 1B; literature findings on LUAD and decreased ADH1B+ CAFs; other finding ADH1B reduces tumor stemness by activating cAMP/PKA/CERB1 signaling.
#KPNA2: Karyopherin alpha 2, nuclear export protein; key target in tumour research, tumour progression, localisation of proteins, poor prognosis
#NFIX: Nuclear factor 1 X-type, member of the nuclear factor I family; transcription factor, expression in adults limited and normally seen in embryonic development, previously found in tumours
#RASGRF1: Ras protein-specific guanine nucleotide-releasing factor 1, RAS-activating guanine exchange factor, fusions common, pathways involving RAS commonly involved in cancer
#TPX2: Targeting protein for Xklp2; microtubule assembly in the M phase, overexpression in various cancers

factor_numeriser <- function(annot) {
  for (i in seq_along(annot)) {
    if (!is.logical(annot[,i])) {
      if (length(unique(annot[,i])) <= 10) {
        annot[,i] <- annot[,i] %>% factor()
      } else {
        out <- tryCatch(
          expr = {temp <- annot[,i] %>% {gsub(",", ".", .)} %>% as.numeric; annot[,i] <- temp},
          warning = function(w) {print(i)},
          error = function(e) {print(i)}
        )
      }
    }
  }
  return(annot)
}

#Correlation Matrix for top 30 DEGs
top_30_cor = exp_df[,which(rank(BH_DGE_df$minimum)<31)] %>% cor(method = "spearman")
pheatmap(top_30_cor, main = "Correlation among the Top 30 differentially expressed genes")

DEG_full_exp = as.data.frame(expression_data) %>% column_to_rownames(colnames(expression_data)[1]) %>% .[expression_data[[1]] %in% unique(top_DEGs),] %>% t() %>% as.data.frame()
Cont_clin_vars = factor_numeriser(anot_df)
selector = vector(length=ncol(Cont_clin_vars))
for (i in seq_along(selector)) selector[i] = is.numeric(Cont_clin_vars[,i])
Cont_clin_vars = Cont_clin_vars[, selector]
Cont_clin_vars = Cont_clin_vars[complete.cases(Cont_clin_vars),]

exp_for_cor <- rownames(DEG_full_exp) %in% rownames(Cont_clin_vars) %>% DEG_full_exp[.,]

cor_genes_clin = function(gene_exp, clin_vars) {
  dat = cbind(gene_exp, clin_vars)
  corr = cor(dat)
  m = ncol(gene_exp)
  n = ncol(Cont_clin_vars)
  corr[1:m, (m + 1):(m + n)]
}

cor_genes_clin(exp_for_cor, Cont_clin_vars) %>% pheatmap()



#DEGs and tumour stage

#We'll use the top 3 DEGs

top_30_DEG = exp_df[,which(rank(BH_DGE_df$minimum)<31)] %>% as.data.frame()
sample_stage = factor_numeriser(anot_df)$ajcc_pathologic_stage
stages = levels(sample_stage) %>% as.vector()

#Prepare results data frame with all the genes and stages to be tested vs. rest
top_30_sample_df = data.frame(matrix(NA, nrow = ncol(top_30_DEG), ncol = n_distinct(stages)))
rownames(top_30_sample_df) = colnames(top_30_DEG)
colnames(top_30_sample_df) = stages


#Run a Wilcox test for each stage vs rest across each gene and save p-value of the test into result dataframe
for (stage in stages){
  not_stage = stages[stages != stage]
  for (gene in colnames(top_30_DEG)) {
    exp_vals = top_30_DEG[[gene]]
    in_group = exp_vals[which(sample_stage == stage)]
    out_group = exp_vals[which(sample_stage %in% not_stage)]
    temporary_test_result = wilcox.test(in_group,out_group,alternative = "two.sided")
    top_30_sample_df[gene,stage] = temporary_test_result$p.value
  }
}

#BH correction
BH_top_30_sample = flatten(top_30_sample_df) %>% p.adjust("BH") %>% matrix(ncol = n_distinct(sample_stage)) %>% as.data.frame(row.names = colnames(top_30_DEG))
colnames(BH_top_30_sample) = stages

pheatmap(BH_top_30_sample, cluster_cols = FALSE, cluster_rows = TRUE, main = "Q-Values of Top 30 genes compared across Subtypes")



