# review posterior parameter estimates

markers <- c("mdr86", "mdr184", "mdr1246", "k13_marcse", "crt76")
models <- c("gneiting_sparse")

posts <- apply(expand.grid(markers, models), 1, function(row){
  out_dir = paste0("output/", paste0(row, collapse = "/"), "/")
  message(out_dir)
  draws <- readRDS(paste0(out_dir, "draws.rds")) %>%
    do.call(what = rbind)
  apply(draws, 2, density)
})

length(posts)
names(posts) <- markers
posts$mdr86

# ugh who knows
df <- data.frame(x = NA, y = NA, pan = NA, mod = NA)
for (i in 1:length(posts$mdr86)){ # 7 - nparams
  for (j in 1:length(posts)){ # models
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



#### prior checks

parameters <- read_rds("output/k13_marcse/gneiting_sparse/parameters.rds")

prior_samps <- lapply(parameters, function(x){
  greta::calculate(x, nsim = 5000)
}) %>%
  as.data.frame()

names(prior_samps) <- c(names(parameters)[1:4], paste("beta", 1:3, sep = "_"))

ggplot(prior_samps %>%
         pivot_longer(cols = everything(), names_to = "param")) +
  geom_density(aes(x = value)) +
  facet_wrap(~param, scales = "free")
  

outdirs <- apply(expand.grid(markers, models), 1, function(row){
  paste0("output/", paste0(row, collapse = "/"), "/")
})

postes <- lapply(outdirs, function(mod){
  message(mod)
  draws <- readRDS(paste0(mod, "draws.rds")) %>%
    do.call(what = rbind)
  mutate(as.data.frame(draws), model = mod)
}) %>%
  do.call(what = rbind)

names(postes) <- c(names(prior_samps), "model")

ggplot(postes %>%
         pivot_longer(cols = -c(model), names_to = "param")) +
  geom_density(aes(x = value, col = model)) +
  facet_wrap(~param, scales = "free")

ggplot(postes %>%
         pivot_longer(cols = -c(model), names_to = "param")) +
  geom_density(aes(x = value, col = model)) +
  geom_density(data = prior_samps %>%
                 pivot_longer(cols = everything(), names_to = "param"),
               aes(x = value)) +
  facet_wrap(~param, scales = "free")










