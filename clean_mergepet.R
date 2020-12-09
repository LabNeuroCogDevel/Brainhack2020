#!/usr/bin/env Rscript
suppressPackageStartupMessages{library(dplry)}

# 20201209WF - init
#   keep only rows we need in giant wide dataframe 

keep_col_regex <- 'lunaid|^age|^visitnum|^sex|vdate|behdatestr|dtbzdate|^fd|uppsp|rist|rt18|tat2|frogET|ysrasr'

read.csv('/Volumes/Phillips/mMR_PETDA/scripts/merged_data.csv') %>%
  select(matches(keep_col_regex)) %>%
  mutate(sesid=glue::glue("{lunaid}_{visitnum}")) %>%
  filter(!is.na(age)) %>%
  write.csv('data/wide.csv', row.names=F, quote=F)

