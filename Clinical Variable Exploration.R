#Clinical Variable Exploration
library(dplyr)
library(tidyverse)
library(ggplot2)

load("/Users/USER/Downloads/LUAD_data.RData")

#Explore pathologic stage data
table(luad_anot_clean$ajcc_pathologic_stage)
ggplot(luad_anot_clean, aes(x=ajcc_pathologic_stage)) + geom_bar() +labs (title="Distribution of pathologic stage", x="pathologic stage", y="frequency")

#Explore Age at diagnosis data
table(luad_anot_clean$paper_Age.at.diagnosis)
#prepare the data
luad_anot_clean$paper_Age.at.diagnosis[luad_anot_clean$paper_Age.at.diagnosis == "[Not Available]"] <-NA
luad_anot_clean$paper_Age.at.diagnosis <-
  as.numeric(as.character(luad_anot_clean$paper_Age.at.diagnosis))

luad_anot_clean %>%
         filter(!is.na(paper_Age.at.diagnosis)) %>% 
         ggplot(aes(x = paper_Age.at.diagnosis)) +
         geom_bar()
#patients are between 41 and 85 years old
mean(luad_anot_clean$paper_Age.at.diagnosis, na.rm = TRUE)
median(luad_anot_clean$paper_Age.at.diagnosis, na.rm = TRUE)

#Explore patient Sex data
table(luad_anot_clean$paper_Sex)
luad_anot_clean$paper_Sex[luad_anot_clean$paper_Sex == "[unknown]"] <-NA

luad_anot_clean %>%
  filter(!is.na(paper_Sex)) %>% 
  ggplot(aes(x = paper_Sex)) +
  geom_bar()
#there are 110 female and 73 male patients

#Explore Smoking status data
table(luad_anot_clean$paper_Smoking.Status)
luad_anot_clean$paper_Smoking.Status[luad_anot_clean$paper_Smoking.Status== "[Not Available]"] <-NA
luad_anot_clean %>%
  filter(!is.na(paper_Smoking.Status)) %>% 
  ggplot(aes(x = paper_Smoking.Status)) +
  geom_bar()

