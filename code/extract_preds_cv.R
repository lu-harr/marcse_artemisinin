# extract predictions for held-out runs
start <- Sys.time()

suppressMessages(source("code/setup.R"))
suppressMessages(source("code/validation_funcs.R"))

args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]
mod <- args[2]
message(paste0("Marker: ", marker))
message(paste0("Model: ", mod))

# here is some run up:
nice_name_lookup <- list("k13_marcse" = "Kelch 13",
                         "crt76" = "Pfcrt-K76T",
                         "mdr86" = "Pfmdr1-N86Y",
                         "mdr184" = "Pfmdr1-Y184F",
                         "mdr1246" = "Pfmdr1-D1246Y")

data_path_lookup <- list("k13_marcse" = "data/clean/moldm_marcse_k13_nomarker.csv",
                         "crt76" = "data/clean/moldm_crt76.csv",
                         "mdr86" = "data/clean/pfmdr_single_mdr86.csv",
                         "mdr184" = "data/clean/pfmdr_single_mdr184.csv",
                         "mdr1246" = "data/clean/pfmdr_single_mdr1246.csv")

folds <- read_rds(paste0("output/", marker, "/", mod, "/cv_folds.rds"))

dat <- extract_preds_cv(data_path = data_path_lookup[[marker]],
                 pred_path = paste0("output/", marker, "/", mod, "/cv_preds/"),
                 folds = folds, 
                 in_buffer = BUFFER,
                 out_buffer = 50000) 
# using large buffer as increased AGG_FACTOR for held-out models
# this completely removed some islands from the raster and chopped up some of the internal 
# waterbodies/holes in the transmission mask
# ultimately we have explicit spatial covariance so it ain't a huge problem if we need to reland 
# some points

write.csv(dat, paste0("output/", marker, "/", mod, "/mut_dat_cv_preds_extracted.csv"),
          row.names = FALSE)

end <- Sys.time()

print(end - start)