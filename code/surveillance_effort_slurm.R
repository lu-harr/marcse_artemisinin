library(terra)
library(sf)
source("code/setup.R")


args <- commandArgs(trailingOnly = TRUE)
snp <- args[1]  # of "k13", "mdr86", "mdr184", "mdr1246", "crt76"
seed <- as.numeric(args[2])
message(paste0("Marker: ", snp))
message(paste0("Seed: ", seed))

out_dir <- paste0(snp, "/")

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

mut_data <- setup_mut_data(in_dat)

message("All years")
wrap_survey_effort(in_dat,
                   paste0(out_dir, "surveillance_effort_", snp, ".grd"),
                   sigma = 1.5, 
                   apply_mask = FALSE)

message("Annual")
wrap_survey_effort(in_dat,
                   paste0(out_dir, "surveillance_effort_", snp, "_annual.grd"),
                   years = unique(mut_data$year),
                   bin_years = FALSE,
                   sigma = 1.5, 
                   apply_mask = FALSE,
                   msg_me = TRUE)
