library(tidyverse)
library(dplyr)

rm(list=ls())

load("LUAD_data.RData")

anot <- luad_anot_clean
anot$paper_patient <- strtrim(anot$barcode, 12)
anot <- anot[, c("barcode", grep("^paper_", colnames(anot), value = TRUE))]

exp <- column_to_rownames(luad_exp_clean, "gene")

# patient_sex <- anot$paper_Sex[!is.na(anot$paper_Sex)]
# names(patient_sex) <- anot$paper_patient[!is.na(anot$paper_Sex)]
# patient_sex <- patient_sex[-which(patient_sex=="[unknown]")]

exp <- as.data.frame(apply(exp, 2, function(x) {
  d <- which(x=="[unknown]" | x=="[not available]" | x=="[Not Available]")
  x[d] <- NA
  return(x)
}))

anot <- as.data.frame(apply(anot, 2, function(x) {
  d <- which(x=="[unknown]" | x=="[not available]" | x=="[Not Available]")
  x[d] <- NA
  return(x)
}))

sum(apply(exp, 2, function(k) sum(is.na(k)))) # no NAs in exp
sum(apply(anot, 2, function(k) sum(is.na(k)))) # 5623 NAs in anot

sum(apply(anot, 2, function(k) sum(is.na(k))))/{dim(anot)[1]*dim(anot)[2]}

length(grep("paper", colnames(anot)))
# How many paper_ columns are retained? >> 20

colnames(anot) <- c(colnames(anot)[1], colnames(anot)[-1] %>% substring(7))
# All data retained apart from staging categorisation

                    