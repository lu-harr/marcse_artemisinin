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

ras <- covariates$pfpr_2000 %>%
  aggregate(fact = 2)

# get mut_data into spatVec
tmp <- mut_data %>%
  vect(geom = c('x', "y"), crs = "epsg:4326") %>%
  terra::project(covariates$pfpr_2000)

ras <- terra::rasterize(tmp, 
                        field = "tested",
                        ras, 
                        fun=sum,
                        background = 0)

gf <- focalMat(ras, 1.5, "Gauss")
test_dens <- raster::focal(ras, gf, pad = TRUE, na.rm=TRUE)
# may take a sec for high res raster ..
# need na.rm=TRUE or we lose edges

# mask
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(test_dens))

test_dens_masked <- mask(test_dens, afr)
plot(test_dens_masked, main="Surveillance effort")
plot(log10(test_dens_masked)) # lol
plot(sqrt(test_dens_masked)) 
# need to have a think about bandwidth, resolution, etc
# have changed resolution: this is now approximately people tested per 100kmsq

# multipanel: time (hist: number tested, number of points)
surveil <- xyFromCell(test_dens_masked, cell = cells(test_dens_masked)) %>%
  as.data.frame() %>%
  mutate(effort = unlist(extract(test_dens_masked, cells(test_dens_masked))))

p1 <- ggplot(data = mut_data %>%
               group_by(year) %>%
               summarise(n = n())) +
  geom_bar(stat = "identity", aes(x = year, y = n)) +
  ylab("Number of locations") +
  xlab("Year") +
  ggtitle("(a)")

p2 <- ggplot(data = mut_data %>%
               group_by(year) %>%
               summarise(tested = sum(tested))) +
  geom_bar(stat = "identity", aes(x = year, y = tested)) +
  ylab("Number of tests") +
  xlab("Year") +
  ggtitle("(b)")

# could just go back to totals? As emphasis is on effort?
# p2 <- ggplot(data = mut_data %>%
#                group_by(year) %>%
#                summarise(present = sum(present),
#                          absent = sum(tested) - present) %>%
#                pivot_longer(!year, names_to = "Tests", values_to = "Count")) +
#   geom_bar(stat = "identity", aes(x = year, y = Count, fill = Tests)) +
#   xlab("Year")

p3 <- ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
  geom_tile(data = surveil, aes(x, y, fill = effort)) +
  scale_fill_viridis_c(na.value = NA, "Tests per \n~100kmsq", trans="sqrt") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("(c)")

library(deeptime)

gg1 <- ggarrange2(p1, p2, layout = rbind(c(1), c(1)), draw = FALSE)
ggarrange2(gg1, p3, widths = c(1,2))
# It only took me 45 mins to work this out I guess





  
