# parallelise model fitting for cross-validation

# with thanks to jerrick:
# https://dept.stat.lsa.umich.edu/~jerrick/courses/stat506_f24/16-parallel-processing.html
start <- Sys.time()

library(parallel)
library(caret)
library(parallelly)
print(paste("Cores:", detectCores()))
print(paste("Cores:", availableCores()))

source("code/setup.R")
source("code/build_design_matrix.R")
source("code/betabinomial_p_rho.R")
source("code/wrap_fit.R")


# commandArgs feels a bit unflashy but I can't keep having separate scripts
args <- commandArgs(trailingOnly = TRUE)
snp <- args[1]  # of "k13", "mdr86", "mdr184", "mdr1246", "crt76"
seed <- as.numeric(args[2])
print(paste0("Marker: ", snp))
print(paste0("Seed: ", seed))

# snp = "k13"
# seed = 123

out_dir <- paste0(snp, "/bb_gne/")

in_dat <- ifelse(snp == "k13",
                 "data/clean/moldm_k13_nomarker.csv",
                 ifelse(snp == "crt76",
                        "data/clean/moldm_crt76.csv",
                        ifelse(snp == "k13_marcse",
                               "data/clean/moldm_marcse_k13_nomarker.csv",
                               paste0(paste0("data/clean/pfmdr_single_", snp, ".csv")))))


print(paste0("Reading in from: ", in_dat))
print("Enforcing min year for surveyor data - 2000")
print("Fitting betabinom")

mut_data <- setup_mut_data(in_dat, min_year = 2000)
write_rds(mut_data, paste0("output/", out_dir, "mut_data.rds"))

NFOLD <- 10
folds <- createFolds(mut_data$present / mut_data$tested, k = NFOLD)
write_rds(folds, paste0("output/", out_dir, "cv_folds.rds"))

# print("para?")
system.time(mclapply(1:NFOLD, function(x){
  fit_betabinom(mut_data = mut_data,
            covariates = covariates,
            pfpr_years = pfpr_years,
            out_dir = out_dir,
            fold = x,
            folds = folds,
            nchains = 6,
            warmup = 3000,
            nsamples = 20000)
}, mc.cores = NFOLD))

# test <- folds[[fold]] # don't need you yet !

