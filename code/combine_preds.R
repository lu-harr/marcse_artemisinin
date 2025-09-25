# combine preds - would be great if this was a function ..
# ugh messing this up - need to make edits at pred step because layer names are wrong

markers <- c("mdr86", "mdr184", "mdr1246", "k13_marcse", "crt76")
models <- c("bb_gne", "gneiting_ahmc")

markers <- c("k13_marcse")
models <- c("gneiting_ahmc")

preds <- apply(expand.grid(markers, models), 1, function(row){
  
  out_dir = paste0("output/", paste0(row, collapse = "/"))#, "/")
  
  # haven't tested this yet
  concat_coverages(out_dir)
  
  to_read <- grep("^20.*\\.grd$", list.files(out_dir), value = TRUE)
  
  tmp = rast(paste0(out_dir, "/", to_read))
  # perhaps trying to write too many things ??
  medians = subset(tmp, grepl("50", names(tmp)))
  sds = subset(tmp, grepl("sd", names(tmp)))
  ci_width = subset(tmp, grepl("97\\.5", names(tmp))) - subset(tmp, grepl("2\\.5", names(tmp)))
  names(ci_width) = paste0(substr(names(ci_width), 1, 5), "CI")
  tmp <- c(medians, sds, ci_width)
  f <- file.path(out_dir, "preds_all.tif")
  terra::writeRaster(tmp, f, overwrite = TRUE, filetype = "GTiff")
  message("done")
})

# mdr86 <- rast("output/mdr184/bb_gne/preds_all.grd")
# plot(mdr86[[c("2006_post_median", "2010_post_median",
#               "2013_post_median")]], breaks = seq(0, 1, 0.01))
# these don't necessarily look like the circ mat outputs ....
# double check model_checks and mock up some figures for Jen
# also re-run fits for circ mat ... I feel like some were falling short

tmp <- rast("output/k13_marcse/bb_gne/preds_all.tif")
looseVis::rast_plot(subset(tmp, grep("50", names(tmp)))) 
looseVis::rast_plot(subset(tmp, grep("CI", names(tmp))))
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/bb_gne/preds_all.tif",
               "k13 gneiting", xlim = c(0, 0.6), ylim = c(0, 0.6))
ggsave("~/Desktop/presentations/marcse/residuals_k13m.png", height = 3, width = 5, scale = 1.5)


tmp <- rast("output/mdr184/bb_gne/preds_all.tif")
rast_plot(subset(tmp, grep("50", names(tmp))))
rast_plot(subset(tmp, grep("CI", names(tmp))))

tmp <- rast("output/mdr1246/bb_gne/preds_all.tif")
rast_plot(subset(tmp, grep("50", names(tmp))))
rast_plot(subset(tmp, grep("CI", names(tmp))))

tmp <- rast("output/crt76/bb_gne/preds_all.tif")
rast_plot(subset(tmp, grep("50", names(tmp))))
rast_plot(subset(tmp, grep("CI", names(tmp))))


# just grab upper_lower to do coverage probs
preds <- apply(expand.grid(markers, models), 1, function(row){
  #message(row)
  out_dir = paste0("output/", paste0(row, collapse = "/"))#, "/")
  #message(out_dir)
  to_read <- grep("^20.*\\.grd$", list.files(out_dir), value = TRUE)
  #message(to_read)
  tmp = rast(paste0(out_dir, "/", to_read))
  # perhaps trying to write too many things ??
  upper = subset(tmp, grepl("97\\.5", names(tmp)))
  lower = subset(tmp, grepl("2\\.5", names(tmp)))
  tmp <- c(upper, lower)
  f <- file.path(out_dir, "preds_upper_lower.tif")
  terra::writeRaster(tmp, f, overwrite = TRUE, filetype = "GTiff")
  message("done")
})



