library(terra)
library(sf)
source("code/setup.R")


args <- commandArgs(trailingOnly = TRUE)
snp <- args[1]  # of "k13", "mdr86", "mdr184", "mdr1246", "crt76"
seed <- as.numeric(args[2])
message(paste0("Marker: ", snp))
message(paste0("Seed: ", seed))

out_dir <- paste0(snp, "/gneiting_ahmc/")

in_dat <- ifelse(snp == "k13",
                 "data/clean/moldm_k13_nomarker.csv",
                 ifelse(snp == "crt76",
                        "data/clean/moldm_crt76.csv",
                        ifelse(snp == "k13_marcse",
                               "data/clean/moldm_marcse_k13_nomarker.csv",
                               paste0(paste0("data/clean/pfmdr_single_", snp, ".csv")))))


message(paste0("Reading in from: ", in_dat))
message(getwd())
message("Enforcing min year for surveyor data - 2000")
mut_data <- setup_mut_data(in_dat, min_year = 2000)
