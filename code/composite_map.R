dat <- read_delim(file = "../moldm/raw/MARCSE_dashboard/TES_outcomes_March_2025.txt", 
                  delim = "\n")
unique(sapply(dat, str_count, ","))

# not sure where the warnings are coming from when I do this:
dat <- read.table(file = "../moldm/raw/MARCSE_dashboard/TES_outcomes_March_2025.txt",
                  sep = ",", header = TRUE) %>%
  mutate(Follow_up = case_when(Follow_up == "Day28" ~ "D28",
                               Follow_up == "Day42" ~ "D42",
                               TRUE ~ Follow_up)) %>%
  rename(x = Long, y = Latitude, year = Start_Year) %>%
  mutate(across(c(PCR_Corrected, year, x, y), ~ as.numeric(.x))) %>%
  mutate(PCR_Corrected = PCR_Corrected/100)

nrow(dat)

ftable(dat$Country)
ftable(dat$Antimalarial)
ftable(dat$Follow_up)
dat$Pf_parasite_positive
dat$Pf_parasite_positive_CI_
dat$PCR_Uncorrected_ACPR
dat$PCR_Corrected # all of them have this ...
dat$Colour_assignment_1


ggplot(dat %>% filter(Antimalarial == "ASAQ")) +
  geom_histogram(aes(x = PCR_Corrected))

ggplot(dat %>% filter(Antimalarial == "AL")) +
  geom_histogram(aes(x = PCR_Corrected))

ggplot(dat %>% filter(Antimalarial == "AL")) +
  geom_point(aes(x = Start_Year, y = PCR_Corrected)) +
  geom_smooth(aes(x = Start_Year, y = PCR_Corrected), method = "lm")

extract_preds <- function(locs,
                           pred_path,
                           ave_tag = "_50",
                           buffer = 0){
  # bring in coords associated with predictions
  preds <- rast(pred_path)
  yrs_pred <- str_extract(names(preds), "\\d{4}")
  
  # get predictions for each row in `locs`
  locs$pred <- NA
  yrs_to_extract <- unique(locs$year)
  for (yr in yrs_to_extract){
    if (yr %in% yrs_pred){
      idx <- which(locs$year == yr)
      val <- terra::extract(preds[[paste0(yr, ave_tag)]], 
                            locs[idx, c("x", "y")],
                            ID = FALSE, search_radius = buffer)
      locs[idx, "pred"] <- val[, paste0(yr, ave_tag)]
    }
  }
  
  locs$pred
}

dat$k13 <- extract_preds(dat, pred_path = "output/k13_marcse/bb_gne/preds_medians.tif")
dat$mdr86 <- extract_preds(dat, pred_path = "output/mdr86/bb_gne/preds_medians.tif")
dat$mdr184 <- extract_preds(dat, pred_path = "output/mdr184/bb_gne/preds_medians.tif")
dat$mdr1246 <- extract_preds(dat, pred_path = "output/mdr1246/bb_gne/preds_medians.tif")
dat$crt76 <- extract_preds(dat, pred_path = "output/crt76/bb_gne/preds_medians.tif")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- cor(unlist(x), unlist(y))
  message(cor(unlist(x), unlist(y)))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

pairs(dat %>% 
        filter(Antimalarial == "AL") %>%
        dplyr::select(PCR_Corrected, year, mdr86, mdr184, mdr1246, crt76, k13) %>%
        rename(start_year = year) %>%
        drop_na() %>%
        as.data.frame(),
      lower.panel = panel.cor
      )

mod <- lm(PCR_Corrected ~ k13 + mdr86 + mdr1246 + mdr184 + crt76, 
          dat %>% filter(Antimalarial == "AL"))
summary(mod)


gmod <- glm(PCR_Corrected ~ k13 + mdr86 + mdr1246 + mdr184 + crt76, 
            data = dat %>% filter(Antimalarial == "AL"),
            family = quasibinomial())
summary(gmod)
# expect -ve coeffs except for 184

tmp <- data.frame(k13 = 0.3, mdr86 = 0.1, mdr184 = 0.55, mdr1246 = 0.1, crt76 = 0.1)
plot(dat$PCR_Corrected, predict(gmod, dat, type = "response"))


library(geosphere)
# Alternatively, could go for recorded prevalences close to TES in time/space?

extract_prevalences <- function(locs, data_path, 
                                # these two params are for spatial and temporal cutoffs:
                                buff = 3.5, year_gap = 3){
  # bring in coords associated with marker prevalences
  mut_data <- setup_mut_data(data_path, min_year = MIN_YEAR) %>%
    st_as_sf(coords = c("x", "y"))
  
  yrs <- unique(locs$year)
  
  out <- data.frame(matrix(NA, nrow = nrow(locs), ncol = 6))
  names(out) <- c("ind", "dist", "present", "tested", "mut_x", "mut_y")
  
  for (yr in yrs){
    mut <- filter(mut_data, year == yr)
    idx <- which(locs$year == yr)
    dists <- st_distance(locs[idx,], mut)
    min_dists <- apply(dists, 1, function(x){
      c(which.min(x), 
        min(x), 
        mut$present[which.min(x)], 
        mut$tested[which.min(x)],
        st_coordinates(mut[which.min(x),])) # points have to be grabbed here if I want them
      }) %>%
      t() %>%
      as.data.frame()
    
    names(min_dists) <- c("ind", "dist", "present", "tested", "mut_x", "mut_y")
    min_dists[min_dists$dist > buff,] <- rep(NA, ncol(min_dists))
    out[idx,] <- min_dists
  }
  
  for (i in 1:year_gap){
    # check adjacent years for unmatched sites
    idx <- which(is.na(out$present))
    if (length(idx) == 0){
      break
    }
    
    dists <- st_distance(locs[idx,], mut_data)
    temps <- outer(locs$year[idx], mut_data$year, function(a, b) abs(a - b))
    for (j in 1:length(idx)){
      # this is a bit gross:
      if (min(dists[j, which(temps[j,] == i)]) < buff){
        cand <- which(temps[j,] == i)[which.min(dists[j, which(temps[j,] == i)])]
        out[idx[j],] <- c(cand, 
                          min(dists[j, which(temps[j,] == i)]), 
                          mut_data$present[cand], 
                          mut_data$tested[cand],
                          st_coordinates(mut_data[cand,]))
      }
    }
  }
  
  out
}


tmp = bind_cols(dat,
                extract_prevalences(st_as_sf(dat, coords = c("x", "y")), 
                    "data/clean/moldm_marcse_k13_nomarker.csv"))

tmp4 <- tmp %>%
  rename(tes_x = x, tes_y = y) %>%
  mutate(idx = 1:nrow(tmp)) %>%
  pivot_longer(cols = matches("_(x|y)$"),
                names_to = c("type", ".value"),
                names_sep = "_")

# have a look at matches:
ggplot() +
  geom_sf(data = afr) +
  geom_line(aes(x = x, y = y, group = idx), data = tmp4) +
  geom_point(aes(x = x, y = y, size = tested, col = present/tested), 
             data = tmp4 %>% filter(type == "mut"), alpha = 0.3) +
  scale_color_viridis_c()

spdat <- st_as_sf(dat, coords = c("x", "y"))
dat <- mutate(dat,
              k13_dat = extract_prevalences(spdat, 
                                            "data/clean/moldm_marcse_k13_nomarker.csv") %>%
                mutate(prev = present/tested) %>%
                dplyr::select(prev) %>%
                unlist(),
              mdr86_dat = extract_prevalences(spdat, 
                                              "data/clean/pfmdr_single_mdr86.csv") %>%
                mutate(prev = present/tested) %>%
                dplyr::select(prev) %>%
                unlist(),
              mdr184_dat = extract_prevalences(spdat, 
                                              "data/clean/pfmdr_single_mdr184.csv") %>%
                mutate(prev = present/tested) %>%
                dplyr::select(prev) %>%
                unlist(),
              mdr1246_dat = extract_prevalences(spdat, 
                                              "data/clean/pfmdr_single_mdr1246.csv") %>%
                mutate(prev = present/tested) %>%
                dplyr::select(prev) %>%
                unlist())
              # ... and crt

pairs(dat %>% 
        filter(Antimalarial == "AL") %>%
        dplyr::select(PCR_Corrected, year, mdr86_dat, mdr184_dat, mdr1246_dat, crt76, k13_dat) %>%
        rename(start_year = year) %>%
        drop_na() %>%
        as.data.frame(),
      lower.panel = panel.cor
)

gmod <- glm(PCR_Corrected ~ mdr86_dat + mdr1246_dat + mdr184_dat + crt76, 
            data = dat %>% filter(Antimalarial == "AL"),
            family = quasibinomial())
summary(gmod)

gmod <- glm(PCR_Corrected ~ k13_dat,
            data = dat, family = quasibinomial())

# this is unfortunately even noisier than the model using map outputs ...

k13 <- rast("output/k13_marcse/bb_gne/preds_medians.tif")
mdr86 <- rast("output/mdr86/bb_gne/preds_medians.tif")
mdr184 <- rast("output/mdr184/bb_gne/preds_medians.tif")
mdr1246 <- rast("output/mdr1246/bb_gne/preds_medians.tif")
crt76 <- rast("output/crt76/bb_gne/preds_medians.tif")

plot(k13$`2024_50` + mdr86$`2024_50` - mdr184$`2024_50` + 
       mdr1246$`2024_50` + crt76$`2024_50`)

plot(k13$`2024_50`)






