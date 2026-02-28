# make up pdps

nice_name_lookup_main <- list("k13_marcse" = "Kelch 13",
                              "crt76" = "Pfcrt-K76T",
                              "mdr86" = "Pfmdr1-N86Y",
                              "mdr184" = "Pfmdr1-Y184F",
                              "mdr1246" = "Pfmdr1-D1246Y")

nice_name_lookup_all <- list("k13_marcse" = "Kelch 13 aggregate",
                             "crt76" = "Pfcrt-K76T",
                             "mdr86" = "Pfmdr1-N86Y",
                             "mdr184" = "Pfmdr1-Y184F",
                             "mdr1246" = "Pfmdr1-D1246Y",
                             "k13snp_A675V" = "Kelch 13 A675V",
                             "k13snp_C469Y" = "Kelch 13 C469Y",
                             "k13snp_P441L" = "Kelch 13 P441L",
                             "k13snp_R561H" = "Kelch 13 R561H",
                             "k13snp_R622I" = "Kelch 13 R622I")

library(tidyverse)

plots = lapply(rep("k13_marcse", 5), function(marker){
  list(pfpr = plot_pdps(out_dir = paste0("output/", marker, "/bb_gne/"), 
                                 target = "pfpr"),
                year = plot_pdps(out_dir = paste0("output/", marker, "/bb_gne/"), 
                                 target = "year"))
})
names(plots) <- names(nice_name_lookup_main)

plot_grid(plots$k13_marcse$pfpr,
      plots$k13_marcse$year,
      plots$crt76$pfpr,
      plots$crt76$year,
      plots$mdr86$pfpr,
      plots$mdr86$year,
      plots$mdr184$pfpr,
      plots$mdr184$year,
      plots$mdr1246$pfpr,
      plots$mdr1246$year, nrow = 5,
      labels = paste0("(", c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), ")"),
      label_fontface = "plain",
      label_y = 1.03, label_x = -0.03, label_size = 12)
ggsave("figures/pdps.png", height = 10, width = 5)
