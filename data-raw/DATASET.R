## code to prepare `DATASET` dataset goes here

library(tidyverse)

# Load the dataset
lrt_yeast <- readr::read_rds("https://github.com/clstacy/GenomicDataAnalysis_Fa23/raw/main/data/yeast_res_edgeR.Rds")

usethis::use_data(lrt_yeast, overwrite = TRUE)


DE_genes <- lrt_yeast %>%
  data.frame() %>%
  filter(PValue < 0.001) %>%
  pull(ORF)

example_enrich_results <- enrichKEGG(gene = DE_genes,
                                     pvalueCutoff = 1,
                                     maxGSSize= Inf,
                                     organism = 'sce'
                                     #options: https://www.genome.jp/kegg/catalog/org_list.html
)

usethis::use_data(example_enrich_results, overwrite = TRUE)

