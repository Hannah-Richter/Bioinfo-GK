#Clinical Variable Exploration
library(dplyr)
library(tidyverse)
library(ggplot2)

load("/Users/USER/Downloads/LUAD_data.RData")

}
table(luad_anot_clean$ajcc_pathologic_stage)
ggplot(luad_anot_clean, aes(x=ajcc_pathologic_stage)) + geom_bar() +labs (title="Distribution of pathologic stage", x="pathologic stage", y="frequency")

table(luad_anot_clean$paper_Age.at.diagnosis)
luad_anot_clean %>%
         filter(!is.na(paper_Age.at.diagnosis)) %>% filter(paper_Age.at.diagnosis !="[Not Available]") %>%
         ggplot(aes(x = paper_Age.at.diagnosis)) +
         geom_bar()