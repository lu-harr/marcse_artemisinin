# review posterior parameter estimates

markers <- c("mdr86", "mdr184", "mdr1246", "k13_marcse", "crt76")
models <- c("gneiting_sparse")

# posts <- apply(expand.grid(markers, models), 1, function(row){
#   out_dir = paste0("output/", paste0(row, collapse = "/"), "/")
#   draws <- readRDS(paste0(out_dir, "draws.rds")) %>%
#     do.call(what = rbind)
#   tmp = lapply(1:ncol(draws), function(x){summary(draws[,x])})
#   names(tmp) = colnames(draws)
#   tmp
# })

posts <- apply(expand.grid(markers, models), 1, function(row){
  out_dir = paste0("output/", paste0(row, collapse = "/"), "/")
  draws <- readRDS(paste0(out_dir, "draws.rds")) %>%
    do.call(what = rbind)
  apply(draws, 2, density)
})

length(posts)
names(posts) <- markers
posts$mdr86

# ugh who knows
df <- data.frame(x = NA, y = NA, pan = NA, mod = NA)
for (i in 1:length(posts$mdr86)){
  for (j in 1:length(posts)){
    tmp <- unlist(posts[[j]][[i]])
    tmp <- data.frame(x = unlist(tmp[grep("x", names(tmp))]),
                      y = unlist(tmp[grep("y", names(tmp))]),
                      pan = names(posts[[j]])[i],
                      mod = names(posts)[j])
    df <- rbind(df, tmp)
  }
}

nrow(df)
df <- df[2:nrow(df),]

ggplot(df) +
  geom_line(aes(x = x, y = y, col = mod)) +
  facet_wrap(~pan, scales = "free")






