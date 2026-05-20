#STEP C Clinical Feature Classification
#define outcome: we want to train our model to predict the expression subtype
unique(luad_anot_clean$paper_expression_subtype)
#TRU=terminal respiratory unit, prox-inflam = proximal-inflammatory, prox-prolif.= proximal-proliferative
#select meaningful features: ajcc_pathologic_stage,paper_sex, paper_Age.at.diagnosis, paper_smoking-status, vlt noch nonsilent mutation
library(dplyr)
library (nnet)

luad_anot_clean_num <- luad_anot_clean %>% 
  select(ajcc_pathologic_stage, paper_Sex, paper_Age.at.diagnosis, paper_Smoking.Status, paper_expression_subtype, paper_Tumor.stage)%>%
  filter(!is.na(paper_Sex))
#convert the categorical values into numerical ones

luad_anot_clean_num$paper_expression_subtype <- factor(luad_anot_clean_num$paper_expression_subtype)
luad_anot_clean_num$paper_Sex <- factor(luad_anot_clean_num$paper_Sex) 
luad_anot_clean_num$paper_Tumor.stage <- factor(luad_anot_clean_num$paper_Tumor.stage)
luad_anot_clean_num$paper_Smoking.Status <- factor(luad_anot_clean_num$paper_Smoking.Status)
luad_anot_clean_num$ajcc_pathologic_stage <- factor(luad_anot_clean_num$ajcc_pathologic_stage)
str(luad_anot_clean_num)

#divide in train and test data 70/30
set.seed(123)
train_d <- createDataPartition(luad_anot_clean_num$paper_expression_subtype, p = 0.7, list = FALSE) 
train_data <- luad_anot_clean_num[train_d, ]
test_data <- luad_anot_clean_num[-train_d, ]

#calculate a multinomal logistic regression (was too complex with all features so I chose 3: Sex, Tumor.stage, Smoking.status)
model_logreg <- multinom(paper_expression_subtype ~ paper_Sex + paper_Tumor.stage + paper_Smoking.Status, data = train_data)
summary(model_logreg)
#choose different features
model_logreg <- multinom(paper_expression_subtype ~ paper_Sex + paper_Age.at.diagnosis + paper_Smoking.Status, data = luad_anot_clean_num)
summary(model_logreg)