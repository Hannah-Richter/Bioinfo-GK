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
         geom_histogram()
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

#Examine Transversion.High.Low 
table(luad_anot_clean$paper_Transversion.High.Low)
#Fusions
table(luad_anot_clean$paper_Fusions)
#EML4-ALK and CLTC-ROS1 are are the most common fusions, CLTC-ROS1 is present in the prox.-inflam expression subtype and EML4-ALK is present in the TRU subtype

#Explore iCluster and CIMP methylation signature via cross-tabulation
table(luad_anot_clean$paper_iCluster.Group, luad_anot_clean$paper_CIMP.methylation.signature.)
#most high CIMp methylation signature samples are in Cluster Group 3 and 4, in Cluster Group 1 are only low CIMP methylation samples
