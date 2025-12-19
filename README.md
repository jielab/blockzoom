### Heritability Estimation and Local LD visualization by blOck

![Fig1](./images/fig1.gif)

![LD](./images/ld-block.png)

source("D:/github/blockzoom/scripts/blockzoom.f.R")
ld_dir <- "D:/data/ldref/csx/ldblk_1kg_eur"
hapref <- "1kg"
glist_file <- NULL
glist <- .read_glist(glist_file)
gwas_file <- "D:/github/blockzoom/test_data/chr1_3blocks.tsv"
block_id <- 31
result <- blockzoom(gwas_file, block_id)