# Corestem project - RMAT Cytokine Analysis.proj
# 20210201 Sehwan Chun at Corestem, Inc.
# 1.1. load environment for analysis

#### 1. Library Loading ####
packs = c("data.table", "readxl", "ggpubr", "writexl", "tidyverse", "caret", "leaps", "MASS")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. Files Loading ####
reg_file = "./Data/Rawdata/Phase12_ALSFRS_12month.txt"
reg_file = read.table(reg_file, header = T)

save.image(file = "./Data/RMAT Cytokine Analysis files.image")
