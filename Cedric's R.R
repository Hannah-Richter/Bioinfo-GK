library(tidyverse)
library(dplyr)

load("LUAD_data.RData")

expression_data = luad_exp_clean


#Histogram to determine logCPM or CPM
expression_data[1,seq(2,length(expression_data[1,]))]%>% unlist(use.names = FALSE) %>% hist(,breaks= 20)

expression_data[,4] %>% unlist(use.names = FALSE) %>% hist(,breaks= 20)



#Histogram not highly skewed -> indicates logCPM over CPM

#Shapiro-Wilks test of all genes separately
norm_test = expression_data[,seq(2,length(expression_data[1,]))] %>% apply(MARGIN =1, function(x) shapiro.test(x))

#Extract p-values from "Shapiro Test" Object and add to list; plot into histogram
Shapiro_p_values = list()
for(i in 1:length(norm_test)) Shapiro_p_values[i]= norm_test[[i]]$p.value
Shapiro_p_values %>% unlist(use.names = FALSE) %>% hist(,breaks=20)
