suppressWarnings(suppressMessages(source("code/setup.R")))
suppressMessages(source("code/build_design_matrix.R"))
source("code/pdp_functions.R")


args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]

message(paste0("Marker: ", marker))

set.seed(0748)

pdp_please(marker)

# doing this in its own script:
# plots <- list(pfpr = plot_pdps(out_dir = paste0("output/", marker, "/bb_gne/"), 
#                                 target = "pfpr"),
#              year = plot_pdps(out_dir = paste0("output/", marker, "/bb_gne/"), 
#                                target = "year"))
# saveRDS(plots, paste0("output/", marker, "/bb_gne/pdps_plotted.rds"))

