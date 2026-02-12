library(tidyverse)

MIN_SAMPLE_SIZE = 10

with_markers <- read.csv("data/clean/moldm_marcse_with_markers.csv") # full

moi <- with_markers %>%
  filter(Marker == "A578S") %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>%
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Present, Site.Name, Country, pubs) %>%
  ungroup() %>%
  filter(Present / Tested <= 1) %>%
  suppressMessages()

ggplot(moi) +
  geom_point(aes(x = year, y = Present/Tested))

wildtypes <- with_markers %>%
  filter(Marker == "wildtype") %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), n=n(), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>% 
  # check how many simultaneous wildtype entries we have?
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country, Present, pubs) %>%
  ungroup() %>%
  suppressMessages()

wildtypes_to_add <- anti_join(
  # locations in `wildtypes` that do not occur in `mutants`,
  # paying attention to `Tested` but NOT to `pubs`
  wildtypes %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country),
  moi %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country)) %>%
  mutate(Present = 0) %>%
  suppressMessages()

with_wildtypes <- full_join(moi, wildtypes_to_add) %>%
  filter(Tested > 5) %>%
  suppressMessages()

to_vis <- with_wildtypes %>% 
  mutate(year_bin = cut(year, breaks = c(min(year) - 1, seq(2009, 2021, 3), max(year)))) %>%
  arrange(Present/Tested, Tested)

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(to_vis, Present == 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, col = "grey50"),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  #new_scale_color() +
  geom_point(data = filter(to_vis, Present > 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_color_manual(name = "Absence", values = c("grey30"), labels=c("")) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt") +
  scale_size_continuous(name = "Tested", range = c(0.2, 5), trans = "sqrt") +
  facet_wrap(~ year_bin, ncol=3) +
  labs(title = "Prevalence of Kelch 13 578S mutations") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_bw() 
ggsave("figures/A578S.png", height = 6, width = 9)


