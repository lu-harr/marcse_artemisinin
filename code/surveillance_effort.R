# here's how we did surveillance effort for JEV
# mosquito_df_expanded <- all_mozzies %>%
#   mutate(n_individuals=case_when(
#     is.na(n_individuals) ~ mean_sample,
#     TRUE ~ round(n_individuals, 0)
#   )) %>% 
#   uncount(n_individuals)
# 
# # mosquito_df_expanded <- data.frame()
# # 
# # # there's a much faster way to do this with uncount
# # # see nat pathogen spatial model script!!
# # for (i_row in 1:nrow(mosquito_df)) {
# #   
# #   # get no. of rows
# #   n_rows <- mosquito_df[i_row,]$n_individuals
# #   
# #   if(is.na(n_rows)){n_rows <- mean_sample}
# #   
# #   info_mat <- mosquito_df[rep(i_row, n_rows),]
# #   
# #   mosquito_df_expanded <- rbind(mosquito_df_expanded, info_mat)
# #   
# # }
# 
# tmp <- rasterize(mosquito_df_expanded[,c('longitude', 'latitude')], ras, background=0, fun='count')
# gf <- focalWeight(ras, mozz_bandwidth_scalar, "Gauss")
# mozz_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)

# so fun would be count tested ... or something along those lines,
# background would be zero

# get mut_data into spatVec
tmp <- mut_data %>%
  vect(geom = c('x', "y"), crs = "epsg:4326") %>%
  terra::project(covariates$pfpr_2000)

ras <- terra::rasterize(tmp, 
                        field = "tested",
                        covariates$pfpr_2000, 
                        fun=sum,
                        background = 0)

gf <- focalWeight(ras, 1.5, "Gauss")
test_dens <- raster::focal(ras, gf, pad = TRUE, na.rm=TRUE)
# may take a sec for high res raster ..
# need na.rm=TRUE or we lose edges

# mask
world_ras <- rasterize(world, covariates$pfpr_2000)
test_dens_masked <- mask(test_dens, world_ras)
plot(test_dens_masked, main="Surveillance effort")
plot(log10(test_dens_masked)) # lol

# need to have a think about bandwidth, resolution, etc




