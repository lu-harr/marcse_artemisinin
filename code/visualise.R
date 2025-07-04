# do some visualisation in here ... 
library(viridisLite)
library(sf)

# mask
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

output_dir <- "output/circmat_pfmdr86/"

###############################################################################
# surveillance effort

survey_effort_panel <- function(in_path, 
                                agg_factor = 1, 
                                pan = "", 
                                main = "",
                                xlab = "Longitude",
                                ylab = "Latitude",
                                lyr_names = main){
  
  test_dens <- lapply(in_path, rast) %>%
    rast()
  if (agg_factor > 1){
    test_dens <- aggregate(test_dens, agg_factor)
  }
  
  # nasty formula to convert to tests per 100kmsq (ish)
  test_dens <- test_dens * 100 / (res(test_dens)[1] * 111) ** 2
  
  # surveil <- xyFromCell(test_dens, cell = cells(test_dens)) %>%
  #   as.data.frame() %>%
  #   cbind(unlist(extract(test_dens, cells(test_dens))))
  surveil <- cbind(xyFromCell(test_dens, cell = cells(test_dens)), 
                   extract(test_dens, cells(test_dens)))
  names(surveil) <- c("x", "y", lyr_names)
  
  surveil <- surveil %>%
    pivot_longer(lyr_names, names_to = "lyr", values_to = "effort") %>%
    mutate(lyr = factor(lyr, levels = lyr_names))
  
  p <- ggplot() +
    geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
    geom_tile(data = surveil, aes(x, y, fill = effort)) +
    geom_sf(data = st_as_sf(afr), colour = "white", fill = NA) +
    scale_fill_viridis_c(na.value = NA, bquote(atop("Tests per","~100"~km^2)), trans="sqrt") +
    facet_wrap(~ lyr) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(paste(pan, main))
  
  p
}

survey_effort_panel(c("output/circmat_k13/surveillance_effort_k13.grd",
                      "output/circmat_crt/surveillance_effort_crt.grd",
                      "output/circmat_pfmdr86/surveillance_effort_pfmdr86.grd"),
                    lyr_names = c("Pfkelch13 (all)", "Pfcrt (76)", "Pfmdr1 (86)"))

# p1 <- survey_effort_panel("output/circmat_k13/surveillance_effort_k13.grd", 
#                           pan = "(a)", main = "Kelch 13", xlab="")
# p2 <- survey_effort_panel("output/circmat_crt/surveillance_effort_crt.grd", 
#                           pan = "(b)", main = "Pfcrt 76", ylab="")
# p3 <- survey_effort_panel("output/pfmdr_hier/surveillance_effort_pfmdr86.grd", 
#                           pan = "(c)", main = "Pfmdr1 86/184/1246", xlab="", ylab="")
# 
# library(patchwork)
# p1 + p2 + p3 +
#   plot_layout(ncol = 3, widths = c(4, 4, 4),
#               axis_titles = "collect")
ggsave("figures/surveillance_effort_all.png", 
        height = 6, width=15, scale=0.7)


# after all that, I might just facet_wrap ....



###############################################################################
# time !




