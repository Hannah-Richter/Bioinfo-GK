library(tidyverse)
library(dplyr)

load("LUAD_data.RData")

luad_anot <- luad_anot_clean
luad_anot$paper_patient <- strtrim(luad_anot$barcode, 12)

luad_exp <- column_to_rownames(luad_exp_clean, "gene")


patient_sex <- luad_anot$paper_Sex[!is.na(luad_anot$paper_Sex)]
names(patient_sex) <- luad_anot$paper_patient[!is.na(luad_anot$paper_Sex)]
patient_sex <- patient_sex[-which(patient_sex=="[unknown]")]

sum(apply(luad_exp, 2, function(k) sum(is.na(k)))) # no NAs in luad_exp

luad_exp <- apply(luad_exp, 2, function(x) {
  d <- which(x=="[unknown]" | x=="[not available]" | x=="[Not Available]")
  x[d] <- NA
  return(x)
})

luad_anot <- apply(luad_anot, 2, function(x) {
  d <- which(x=="[unknown]" | x=="[not available]" | x=="[Not Available]")
  x[d] <- NA
  return(x)
})
