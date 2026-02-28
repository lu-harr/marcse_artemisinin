#' Make predictions to arbitrary coords
#'
#' @param coords data.frame with columns in `coord_cols` and `design_cols`
#' @param draws greta_mcmc_list
#' @param parameters named list of greta arrays
#' @param random_field greta array
#' @param nsim numeric
#' @param coord_cols vector of chars
#' @param design_cols vector of chars
#' @param out chars: location to write to
#'
#' @returns 0
#' @export
#'
#' @examples
predict_for_pdp <- function(coords, 
                            draws, 
                            parameters,
                            random_field,
                            nsim = 100,
                            coord_cols = c("x", "y", "year"),
                            design_cols = c("year"),
                            out = ""){
  
  # finish off design matrix
  X_pred <- coords %>%
    drop_na() %>%
    mutate(intercept = 1)
  
  # project random field to coordinates we would like predictions for
  random_field_pred <- greta.gp::project(random_field, X_pred[,coord_cols])
  
  mut_freq <- (X_pred[,design_cols] %*% parameters$beta + 
                 random_field_pred) %>%
    ilogit()
  
  post_sims <- greta::calculate(mut_freq,
                                values = draws,
                                nsim = nsim,
                                trace_batch_size = 1) # reducing: will take longer, use less mem
  
  message("predicted")
  write_rds(post_sims, out)
}


#' Make up sims for partial dependence plots of PfPR and year
#'
#' @param marker character, e.g. "k13_marcse"
#'
#' @returns 0
#' @export
#'
#' @examples
pdp_please <- function(marker){
  # for given marker, please give me predictions along sets of time- and pfpr-varying coords
  out_dir <- paste0("output/", marker, "/bb_gne/")
  
  draws <- readRDS(paste0(out_dir, "draws.rds"))
  random_field <- readRDS(paste0(out_dir, "random_field.rds"))
  parameters <- readRDS(paste0(out_dir, "parameters.rds"))
  
  pfpr_range <- minmax(pfpr) %>% range()
  pfpr_to_show <- seq(pfpr_range[1], pfpr_range[2], length.out = 100)  
  years_to_show <- 2000:2028
  
  mut_data <- setup_mut_data(data_path_lookup[[marker]], 
                             min_year = MIN_YEAR, buffer = 5000)
  tmp <- build_design_matrix(covariates = pfpr,
                             degs_to_rads = TRUE,
                             coords = mut_data[,c("x", "y", "year")],
                             temporal_var = TRUE,
                             temporal_covt_range = 2000:2028)
  scaled_years <- tmp$scaled_years
  scaled_years <- extended_scaled_years(scaled_years, c(2000, 2025:2028))
  
  mut_data <- dplyr::select(tmp$df, x, y, x_rd, y_rd, year, year_scaled, pfpr) %>%
    drop_na() # some NA pfprs
  
  # make up an X_pred, ~ snag here: what values do I use for other columns?
  # make predictions, given draws from posterior
  # show on one axis ...
  coords_pfpr <- mut_data %>%
    dplyr::select(-c(pfpr)) %>%
    merge(data.frame(pfpr = pfpr_to_show))
  
  write.csv(coords_pfpr, paste0(out_dir, "coords_pfpr.csv"), row.names = FALSE)
  
  predict_for_pdp(coords_pfpr,
                  draws,
                  parameters,
                  random_field,
                  coord_cols = c("x_rd", "y_rd", "year_scaled"),
                  design_cols = c("intercept", "year_scaled", "pfpr"),
                  nsim = 500,
                  out = paste0(out_dir, "pdp_pfpr.rds"))
  
  coords_year <- mut_data %>%
    dplyr::select(-c(year, year_scaled)) %>%
    merge(data.frame(year = names(scaled_years),
                     year_scaled = unlist(scaled_years)))
  
  write.csv(coords_year, paste0(out_dir, "coords_year.csv"), row.names = FALSE)
  
  predict_for_pdp(coords_year,
                  draws,
                  parameters,
                  random_field,
                  coord_cols = c("x_rd", "y_rd", "year_scaled"),
                  design_cols = c("intercept", "year_scaled", "pfpr"),
                  nsim = 500,
                  out = paste0(out_dir, "pdp_year.rds"))
  
  message("all done")
  return(0)
}


#' Make up ggplot partial dependence plots given outputs written by pdp_please()
#'
#' @param out_dir character. directory to read from and write to
#' @param target character. "pfpr" or "year"
#'
#' @returns
#' @export
#'
#' @examples
plot_pdps <- function(out_dir, target){
  pdps <- read_rds(paste0(out_dir, "pdp_", target, ".rds"))
  
  meds <- apply(pdps$mut_freq[,,1], 2, median)
  
  coords <- read.csv(paste0(out_dir, "coords_", target, ".csv")) %>%
    mutate(pred = meds,
           id = paste(x,y,year))
  
  coords_quants <- coords %>%
    group_by(!! sym(target))  %>%
    summarise(q = list(quantile(pred, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
    unnest_wider(q) %>%
    ungroup()
  
  ggplot(data = coords) +
    geom_line(aes(group = id, x = !! sym(target), y = pred), linewidth = 0.1, col = "grey50") +
    geom_ribbon(data = coords_quants,
                aes(x = pfpr, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.4, fill = "#9AD3EC") +
    geom_line(data = coords_quants,
              aes(x = pfpr, y = `50%`),
              col = "#0986B3", linewidth = 1.5) +
    xlab(ifelse(target == "pfpr", "PfPR (scaled)", "Year")) +
    ylab("Prevalence") +
  theme_bw()
}






