#STEP C Clinical Feature Classification
#define outcome: we want to train our model to predict the expression subtype
unique(luad_anot_clean$paper_expression_subtype)
#TRU=terminal respiratory unit, prox-inflam = proximal-inflammatory, prox-prolif.= proximal-proliferative
#select meaningful features: ajcc_pathologic_stage,paper_sex, paper_Age.at.diagnosis, paper_smoking-status, vlt noch nonsilent mutation
library(dplyr)
library (nnet)
#remove rows that are only NAs
luad_anot_clean_num <- luad_anot_clean %>% 
  select(ajcc_pathologic_stage, paper_Sex, paper_Age.at.diagnosis, paper_Smoking.Status, paper_expression_subtype, paper_Tumor.stage)%>%
  filter(!is.na(paper_Sex))
#look how many single NAs are left, there are 19 NA values 
colSums(is.na(luad_anot_clean_num))
luad_anot_clean_num<- na.omit(luad_anot_clean_num)

#prepare expression_subtype data, convert names to valid R variable names
luad_anot_clean_num$paper_expression_subtype[luad_anot_clean_num$paper_expression_subtype == "prox.-inflam"] <- "PI"
luad_anot_clean_num$paper_expression_subtype[luad_anot_clean_num$paper_expression_subtype == "prox.-prolif."] <- "PP"
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
model_logreg1 <- multinom(paper_expression_subtype ~ paper_Sex + paper_Tumor.stage + paper_Smoking.Status, data = train_data)
summary(model_logreg)
#choose different features
model_logreg2 <- multinom(paper_expression_subtype ~ paper_Sex + paper_Age.at.diagnosis + paper_Smoking.Status, data = luad_anot_clean_num)
summary(model_logreg)
#this model with these 3 features makes the most sense to predict expression_subtype (model has the highest accuracy)
model_logreg <-multinom(paper_expression_subtype ~ paper_Tumor.stage + paper_Age.at.diagnosis + paper_Smoking.Status, data = luad_anot_clean_num)
#summary(model_logreg)

#predicting with the test data
pred_logreg <- predict(model_logreg, newdata = test_data)
library(caret)

confusionMatrix(
  data = pred_logreg,
  reference = test_data$paper_expression_subtype)

#random forest
library(randomForest)
model_randomforest <- randomForest(paper_expression_subtype ~ ., data = train_data, importance = TRUE, ntree = 500)
#check for the most important features 
varImpPlot(model_randomforest, main = "Feature Importance")
#predict with the test data
pred_randomforest <- predict(model_randomforest, test_data)
confusionMatrix(pred_randomforest, test_data$paper_expression_subtype)

#k-nearest neighbors
#first normalize the data (only features not the outcome, expression_subtype is in column 5 so we are excluding it)
data_preproc <- preProcess(train_data[, -5], method = c("center", "scale"))
train_norm <- predict(data_preproc, train_data[, -5])
test_norm <- predict(data_preproc, test_data[, -5])

train_norm$paper_expression_subtype<-train_data$paper_expression_subtype
test_norm$paper_expression_subtype<-test_data$paper_expression_subtype

#train the model
set.seed(123)
model_knn <- train(
  paper_expression_subtype ~ .,
  data = train_norm,
  method = "knn",
  tuneGrid = expand.grid(k = 1:10),
  trControl = trainControl(method = "cv", number = 5, classProbs = TRUE)
)
#testing the knn model
pred_knn <- predict(model_knn, newdata = test_norm)
confusionMatrix(pred_knn, test_norm$paper_expression_subtype)