library(tidyverse)
library(dplyr)
library(randomForest)
library(ggplot2)

rm(list=ls())

load("LUAD_data.RData")

NA_flattener <- function(data) {
  as.data.frame(apply(data, 2, function(x) {
    x[which(x=="[unknown]" | x=="[not available]" | x=="[Not Available]")] <- NA
    return(x)
  }))
}

NA_counter <- function(data, per_col=F) {
  s <- apply(data, 2, function(k) sum(is.na(k)))
  if (per_col) return(s) else return (sum(s))
}

level_counter <- function (data) {
  return(sapply(levels(data), function(x) sum(data==x)))
}

anot <- luad_anot_clean
anot$paper_patient <- strtrim(anot$barcode, 12)
anot <- anot[, c("barcode", grep("^paper_", colnames(anot), value = TRUE))]

exp <- column_to_rownames(luad_exp_clean, "gene")

# patient_sex <- anot$paper_Sex[!is.na(anot$paper_Sex)]
# names(patient_sex) <- anot$paper_patient[!is.na(anot$paper_Sex)]
# patient_sex <- patient_sex[-which(patient_sex=="[unknown]")]

exp <- NA_flattener(exp)

anot <- NA_flattener(anot)

length(grep("paper", colnames(anot)))
# How many paper_ columns are retained? >> 20

colnames(anot) <- c(colnames(anot)[1], colnames(anot)[-1] %>% substring(7))
# All data retained apart from staging categorisation

NA_counter(exp) # no NAs in exp
NA_counter(anot) # 5623 total NAs in anot
NA_counter(anot, T) # looking at distribution across rows

anot %>% NA_counter(T) %>% .['barcode']
# Just checking the 

anot$expression_subtype <- anot$expression_subtype %>% factor()
anot <- anot %>% mutate(is.tumour = substring(barcode, 14, 15)=="01") 
# based on the number of barcode segment 4, encoding tumour as 01 and healty as 11

train.anot <- anot[!is.na(anot$expression_subtype),]
train.exp <- exp[,match(train.anot$barcode, colnames(exp))]

exp.var <- train.exp %>% apply(1, var)
t50.train.exp.var <- t(train.exp[names(sort(exp.var)[1:50]),])
train <- data.frame(factor(train.anot$expression_subtype), t50.train.exp.var)
colnames(train) <- c("expression_subtype", colnames(train)[-1])

rf.reg <- randomForest(expression_subtype ~ ., data = train, importance = TRUE, ntree = 500)

# log.reg <- glm(train$expression_subtype ~ t50.train.exp.var, family="binomial")

predict.anot <- anot[is.na(anot$expression_subtype),]
predict.exp <- exp[,match(predict.anot$barcode, colnames(exp))]

t50.predict.exp.var <- t(predict.exp[names(sort(exp.var)[1:50]),])
predict <- data.frame(factor(predict.anot$expression_subtype), t50.predict.exp.var)
colnames(predict) <- c("expression_subtype", colnames(predict)[-1])


rf.pred <- predict(rf.reg, predict)

# log.pred <- predict(log.reg, predict)

ggplot() +
  geom_bar(aes(x=rf.pred), color='black', fill='white') +
  labs(title='Distribution in Random Forest Prediction') +
  theme_bw()

ggplot() +
  geom_bar(aes(x=train$expression_subtype), color='black', fill='white') +
  labs(title='Distribution in Training Data') +
  theme_bw()

level_counter(train$expression_subtype) - level_counter(rf.pred)