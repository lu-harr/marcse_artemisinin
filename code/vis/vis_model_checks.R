# vis model checks
# (trace plots, marginal posterior densities)

library(GGally)
library(brms)
library(bayesplot)

output_dir <- "output/mdr1246/circmat/"
out_tag <- "k13_gneiting_sparse"

draws <- read_rds(paste0(output_dir, "draws.rds"))

r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)

# post <- as_draws_df(draws) %>%
#   rename("Lengthscale (spatial)" = "circmat_len",
#          "Variance (spatial)" = "circmat_var",
#          "Lengthscale (temporal)" = "expo_len",
#          "Variance (temporal)" = "expo_var", # (years are scaled - unscale relevant params?)
#          "Nugget variance" = "nugget_var",
#          "Beta (intercept)" = "beta[1,1]",
#          "Beta (scaled year)" = "beta[2,1]",
#          "Beta (PfPR)" = "beta[3,1]")

post <- as_draws_df(draws) %>%
  rename("Lengthscale (spatial)" = "gneiting_len",
         "Std Deviation (spatiotemporal)" = "gneiting_sd",
         "Lengthscale (temporal)" = "gneiting_tim",
         "Nugget std dev" = "nugget_sd",
         "Beta (intercept)" = "beta[1,1]",
         "Beta (scaled year)" = "beta[2,1]",
         "Beta (PfPR)" = "beta[3,1]")

color_scheme_set("purple") # the bayesplot scheme is much better suited to this application ..
bayesplot::mcmc_trace(post)
ggsave(paste0("figures/", out_tag, "_trace.png"), height=10, width=15)

modified_density = function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + 
    scale_fill_manual(values = c("#e5cce5","#bf7fbf","#a64ca6","#800080","#660066","#400040"))
}

# ggplot is the silliest darn software going
# why would anyone want or need to change a colour palette?
# beats me
modified_points = function(data, mapping, ...) {
  ggally_points(data, mapping, ...) + 
    scale_color_manual(values = c("#e5cce5","#bf7fbf","#a64ca6","#800080","#660066","#400040"))
}

modified_cor = function(data, mapping, ...){
  ggally_cor(data, mapping, ...) +
    scale_color_manual(values = c("#e5cce5","#bf7fbf","#a64ca6","#800080","#660066","#400040"))
}

# love the purps but they're a bit too light to use for cor
tmp <- post %>%
  mutate(chain = as.factor(.chain)) %>%
  ggpairs(columns = 1:6,
          mapping = aes(colour = chain),
          lower = list(continuous = wrap(modified_points, alpha = 0.4)), # can't seem to turn alpha down here?
          diag = list(continuous = wrap(modified_density, alpha = 0.5)),
          upper = list(continuous = modified_cor)) # probably don't need corrs for chains?
ggsave(paste0("figures/chain_corr_",out_tag,".png"), plot = tmp, height=12, width=12)

