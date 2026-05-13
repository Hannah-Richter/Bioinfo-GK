library(tidyverse)
library(dplyr)

load("LUAD_data.RData")

expression_data = luad_exp_clean

#row1 = expression_data %>% dplyr::slice(1)
#values = row1 %>% dplyr::select(-gene) %>% unlist(use.names = FALSE)


#Histogram to determine logCPM or CPM
expression_data[1,seq(2,length(expression_data[1,]))]%>% unlist(use.names = FALSE) %>% hist(,breaks= 20)


#Histogram not highly skewed -> indicates logCPM over CPM

  
norm_test = expression_data[,seq(2,length(expression_data[1,]))] %>% apply(MARGIN =1, function(x) shapiro.test(x))

Shapiro_p_values = c()
for(i in length(norm_test)) Shapiro_p_values[i] = norm_test[[1]]$p.value

                                        
                                        
                                        
matri = norm_test %>% unlist() %>% matrix(., nrow=4)
hist(matri[1,])
