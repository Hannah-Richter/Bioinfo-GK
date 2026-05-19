library(tidyverse)
library(dplyr)
library(randomForest)
library(ggplot2)
library(glmnet)

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

in_out_transformer <- function (vector, inner) {
  if ((!inner %in% unique(vector))) stop("There's no position in your vector equivalent to the inner you provided.")
  return (factor(sapply(vector, function(x) {
    if (x == inner) return(inner) else return (paste0("not.", inner))
  })))
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

anot %>% NA_counter(T) %>% .["barcode"]
# Just checking the NAs in barcode column

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

# Binomial not possible with three levels. Hence splitting into in- and out-group
for (i in levels(train$expression_subtype)) {
  train.in <- train
  train.in$expression_subtype <- in_out_transformer(train.in$expression_subtype, i)
  assign(paste0("log.reg.", i), glm(expression_subtype ~ ., data = train.in, family="binomial"))
}



prediction.anot <- anot[is.na(anot$expression_subtype),]
prediction.exp <- exp[,match(prediction.anot$barcode, colnames(exp))]

t50.predict.exp.var <- t(prediction.exp[names(sort(exp.var)[1:50]),])
prediction <- data.frame(factor(prediction.anot$expression_subtype), t50.predict.exp.var)
colnames(prediction) <- c("expression_subtype", colnames(prediction)[-1])


rf.pred <- predict(rf.reg, prediction)

for (i in levels(train$expression_subtype)) {
  assign(paste0("log.pred.", i), factor(ifelse(predict(get(paste0("log.reg.", i)), prediction) > 0.5, i, paste0("not.", i))))
}

ggplot() +
  geom_bar(aes(x=rf.pred), color="black", fill="white") +
  labs(title="Distribution in Random Forest Prediction") +
  theme_bw()

ggplot() +
  geom_bar(aes(x=train$expression_subtype), color="black", fill="white") +
  labs(title="Distribution in Training Data") +
  theme_bw()

for (i in levels(train$expression_subtype)) {
  show(ggplot() +
    geom_bar(aes(x=get(paste0("log.pred.", i))), color="black", fill="white") +
    labs(title=paste("Distribution in Logistic Prediction with", i, "vs.", paste0("not-", i)), x=paste0("log.pred.", i)) +
    theme_bw())
}

level_counter(train$expression_subtype) - level_counter(rf.pred)
