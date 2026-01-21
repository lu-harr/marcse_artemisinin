# script for validation outputs

# for nearest neighbour index - from Foo and Flegg
library(tensorflow)
source("code/setup.R")
source("code/build_design_matrix.R") # for year scaling

# bringing in some backend functions from greta.gp:
source("~/greta.gp.st.on.earth/R/tf_kernels.R")
source("code/betabinomial_p_rho.R")

library(looseVis)
library(iddoPal)
library(cowplot)

source("code/validation_funcs.R")

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

bb_paths <- lapply(names(nice_name_lookup),
                   function(marker){paste0("output/", marker, "/bb_gne/")})
names(bb_paths) <- names(nice_name_lookup)

# read in mut_data and associate each record with predicted prevalence for 
# relevant year:
mut_dat_assoc_with_preds <- lapply(names(nice_name_lookup), function(marker){
  extract_preds(data_path = data_path_lookup[[marker]],
                pred_path = paste0(bb_paths[[marker]], "preds_medians.tif"))
}) %>%
  setNames(names(nice_name_lookup)) %>%
  suppressMessages()

mut_dat_assoc_with_preds_cv <- lapply(names(nice_name_lookup), function(marker){
  read.csv(paste0(bb_paths[[marker]], "mut_dat_cv_preds_extracted.csv"))
}) %>%
  setNames(names(nice_name_lookup))

rmses <- lapply(mut_dat_assoc_with_preds, function(x){rmse(x)})
rsq <- lapply(mut_dat_assoc_with_preds, function(x){unadjusted_rsq(x)})


# given predicted prevalence, take posterior samples at location of all observations
# and compare quantiles of samples to observed number of cases with marker
sim_coverages <- lapply(names(nice_name_lookup), function(marker){
  message(marker)
  coverage_probabilities_from_observation_model(mut_dat_assoc_with_preds[[marker]],
                                                bb_paths[[marker]],
                                                probs = seq(0, 1, 0.01),
                                                nsim = 500)
})
# (all of those samples and their summarisation takes a bit of a while)
names(sim_coverages) <- names(nice_name_lookup)

sim_coverages <- lapply(names(sim_coverages), function(x){
  mutate(sim_coverages[[x]], marker = x)
}) %>%
  do.call(what = rbind) %>%
  mutate(recs = "All records") %>%
  bind_rows(
    coverage_probabilities_from_observation_model(
      mut_dat_assoc_with_preds[["k13_marcse"]] %>%
        filter(present > 0),
      bb_paths[[marker]],
      probs = seq(0, 1, 0.01),
      nsim = 500) %>%
                mutate(recs = "Presences only",
                       marker = "k13_marcse")
    )



posterior_predictive_ecdfs <- lapply(names(nice_name_lookup), function(marker){
  message(marker)
  posterior_predictive_check(mut_dat_assoc_with_preds[[marker]],
                             bb_paths[[marker]])
})
names(posterior_predictive_ecdfs) <- names(nice_name_lookup)
posterior_predictive_ecdfs <- lapply(names(posterior_predictive_ecdfs), function(x){
  mutate(posterior_predictive_ecdfs[[x]], marker = x)
}) %>%
  do.call(what = rbind)


p1 <- coverages_fig(paste0("output/", c("k13_marcse", "crt76", "mdr86", 
                                        "mdr184", "mdr1246"), "/bb_gne/"))

p2 <- ggplot(data = sim_coverages %>%
               mutate(marker = factor(unlist(nice_name_lookup[marker]),
                                      levels = nice_name_lookup)), 
             aes(x = widths, y = cover)) +
  # geom_point(aes(group = interaction(marker, recs), col = marker), size = 0.5) +
  geom_line(aes(group = interaction(marker, recs), col = marker, linetype = recs)) +
  geom_abline(slope = 1, col = "grey") +
  scale_color_manual(values = viridis(5)) +
  xlab("Posterior predictive interval width") +
  ylab("Coverage probability") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  theme_bw()

p3 <- ggplot(bind_rows(posterior_predictive_ecdfs %>% 
                         filter(name == "leq") %>%
                         mutate(dat = "All records"), 
                       posterior_predictive_ecdfs %>% 
                         filter(name == "leq" & marker == "k13_marcse" & present != 0) %>%
                         mutate(dat = "Presences only")) %>%
               mutate(marker = factor(unlist(nice_name_lookup[marker]),
                                      levels = nice_name_lookup)),
             aes(x = value)) +
  stat_ecdf(geom = "step", 
            aes(group = interaction(marker, dat), 
                col = marker, linetype = dat)) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  scale_color_manual("Marker", values = viridis(5)) +
  scale_linetype("") +
  ylab("ECDF") +
  xlab("Pr(posterior predictive samples <= observed data)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()
# Kelch 13: high number of zeroes and low probabilities: 100% of samples leq observation

leg = get_legend(p3)
plot_grid(
  plot_grid(p1 + theme(legend.position = "none"),
            p2 + theme(legend.position = "none"), 
            p3 + theme(legend.position = "none"), ncol = 1, align = "v",
            labels = c("(a)", "(b)", "(c)"), label_fontface = "plain", label_x = 0.11, label_y = 0.98),
  leg, nrow = 1, rel_widths = c(1, 0.3)
)
ggsave("figures/coverages.png", width = 7, height = 7)






# mut_data <- mut_data %>%
#   mutate(nn = nn_measure(mut_data, 
#                          draws_path),
#          nnplot = 2**nn,
#          nnplot = nnplot/max(nnplot))
# 
# ggplot(mut_data) +
#   geom_point(aes(x = x, y = y, col = nnplot))
# 
# ggplot(mut_data) +
#   geom_point(aes(x = present/tested, y = abs(pred - present/tested), col = nnplot)) +
#   scale_color_viridis_c("Mean distance\nto other points")
# 
# ggplot(mut_data) +
#   geom_sf(data = afr) +
#   geom_point(aes(x = x, y = y, col = nnplot), alpha = 0.3) +
#   scale_color_viridis_c("Mean distance\nto other points") +
#   xlab("Longitude") +
#   ylab("Latitude")

# # e.g.:
# # might want to re-land some points inside of model fitting
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/bb_gne/preds_medians.tif",
               xlim = c(0, 0.6), ylim = c(0, 0.6),
               ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
# tmp = extract_preds("data/clean/moldm_marcse_k13_nomarker.csv",
#               "output/k13_marcse/bb_gne/preds_medians.tif", ave_tag = "_50", 
#               buffer = 100000)

# obs_prev_panel_nn("data/clean/moldm_marcse_k13_nomarker.csv",
#                   "output/k13_marcse/gneiting_sparse/preds_medians.tif",
#                   "output/k13_marcse/gneiting_sparse/",
#                   xlim = c(0, 0.6), ylim = c(0, 0.6),
#                   buffer = 100000)

tmp <- lapply(names(data_path_lookup), function(marker){
  message(marker)
  obs_prev_panel(data_path_lookup[[marker]],
                 paste0(bb_paths[[marker]], "preds_medians.tif"),
                 xlim = c(0, 1), ylim = c(0, 1),
                 ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
})
names(tmp) <- names(data_path_lookup)

# grid these into supp:
p <- plot_grid(tmp$k13_marcse, tmp$crt76, tmp$mdr86, tmp$mdr184, tmp$mdr1246,
          nrow = 3, ncol = 2)
ggsave("figures/obs_prev_panelled.png", p, height = 9, width = 8, scale = 1.4)
# legends will need a fiddle
# could be 3 * 5 with UGA inset ...








