start <- Sys.time()

suppressMessages(source("code/setup.R"))
suppressMessages(source("code/predict_to_raster.R"))

args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]
mod <- args[2]
message(paste0("Marker: ", marker))
message(paste0("Model: ", mod))

# set this to the location where all the inference outputs are at:
# out_dir <- "output/mdr86/gneiting_sparse/"
out_dir <- paste0("output/", marker, "/", mod, "/")

concat_coverages(out_dir)
concat_preds(out_dir, 
             medians = TRUE, 
             sds = TRUE, 
             sdscaled = TRUE, 
             ciwidth = TRUE,
             upper = TRUE,
             lower = TRUE)

end <- Sys.time()

print(end - start)