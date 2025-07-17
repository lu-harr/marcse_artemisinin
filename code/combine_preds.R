# combine preds

markers <- c("mdr86", "mdr184", "mdr1246")
models <- c("gneiting_sparse")

preds <- apply(expand.grid(markers, models), 1, function(row){
  out_dir = paste0("output/", paste0(row, collapse = "/"), "/")
  #message(out_dir)
  to_read <- grep("^20.*\\.grd$", list.files(out_dir), value = TRUE)
  #message(to_read)
  writeRaster(x = rast(paste0(out_dir, to_read)), 
              filename = paste0(out_dir, "preds_all.grd"), overwrite = TRUE)
})

mdr86 <- rast("output/mdr184/gneiting_sparse/preds_all.grd")
plot(mdr86[[c("2006_post_median", "2010_post_median",
              "2013_post_median")]], breaks = seq(0, 1, 0.01))
# these don't necessarily look like the circ mat outputs ....
# double check model_checks and mock up some figures for Jen
# also re-run fits for circ mat ... I feel like some were falling short