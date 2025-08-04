# combine preds - would be great if this was a function ..
# ugh messing this up - need to make edits at pred step because layer names are wrong

markers <- c("mdr86", "mdr184", "mdr1246", "k13", "crt76")
markers <- c("mdr184", "mdr1246", "crt76")
models <- c("gneiting_sparse")

preds <- apply(expand.grid(markers, models), 1, function(row){
  out_dir = paste0("output/", paste0(row, collapse = "/"), "/")
  #message(out_dir)
  to_read <- grep("^20.*\\.grd$", list.files(out_dir), value = TRUE)
  message(to_read)
  tmp = rast(paste0(out_dir, to_read))
  message(names(tmp))
  # get rid of this !
  names(tmp) = paste0(rep(2000:2024, each = 2), c("_post_median", "_post_sd"))
  writeRaster(x = tmp, 
              filename = paste0(out_dir, "preds_all.grd"), overwrite = TRUE)
})

mdr86 <- rast("output/mdr184/gneiting_sparse/preds_all.grd")
plot(mdr86[[c("2006_post_median", "2010_post_median",
              "2013_post_median")]], breaks = seq(0, 1, 0.01))
# these don't necessarily look like the circ mat outputs ....
# double check model_checks and mock up some figures for Jen
# also re-run fits for circ mat ... I feel like some were falling short