# need to grab merged MARCSE data here ....
# not sure what KB has in here

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
library(terra)
library(iddoPal)
library(patchwork)

marker_reference <- readxl::read_xlsx("data/marker_index.xlsx") %>%
  rbind(c("A578A/S", "of interest"),
        c("A578S", "of interest"))

YEAR_LOWER_BOUND <- 2000
MIN_SAMPLE_SIZE <- 10

theme_set(theme_bw())

world <- ne_countries(scale="medium", returnclass = "sf")
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

# moldm <- read_xlsx("raw/SO_attachment_20250610/Surveyor data 250609.xlsx") %>%
#   filter(!is.na(`Start Year`) & !is.na(`End Year`)) %>% # remove where both are NA
#   filter(`End Year` > 1960 & `End Year` < 2500) %>%
#   filter(`Start Year` > 1960 & `Start Year` < 2500) %>%
#   filter(grepl("k13", Marker)) %>%
#   mutate(strip_marker = gsub("k13 ", "", Marker)) %>% # strip k13
#   left_join(marker_reference, 
#             by = join_by(strip_marker == marker)) %>%
#   mutate(year = round((`Start Year` + `End Year`) / 2, 0),
#          # if one or the other is not complete, populate with the value we have:
#          year = case_when(is.na(year) & !is.na(`Start Year`) ~ `Start Year`,
#                           is.na(year) & !is.na(`End Year`) ~ `End Year`,
#                           TRUE ~ year),
#          mutant = !is.na(status) & status != "Not associated") %>%
#   filter(Continent == "Africa") %>%
#   suppressWarnings()

raw_moldm <- function(path){
  read.csv(path) %>%
    mutate(across(c(Start.Year, End.Year, Present, Tested), as.numeric)) %>% 
    filter(!is.na(Start.Year) & !is.na(End.Year)) %>% # remove where both are NA
    filter(End.Year > 1960 & End.Year < 2500) %>%
    filter(Start.Year > 1960 & Start.Year < 2500) %>%
    filter(grepl("[k|K]13", Marker)) %>%
    mutate(strip_marker = gsub("[k|K]13 ", "", Marker)) %>% # strip k13
    left_join(marker_reference, 
              by = join_by(strip_marker == marker)) %>%
    mutate(year = round((Start.Year + End.Year) / 2, 0),
           # if one or the other is not complete, populate with the value we have:
           year = case_when(is.na(year) & !is.na(Start.Year) ~ Start.Year,
                            is.na(year) & !is.na(End.Year) ~ End.Year,
                            TRUE ~ year),
           mutant = !is.na(status) & status != "Not associated" & status != "of interest") %>%
    filter(Continent == "Africa") %>%
    filter(year >= YEAR_LOWER_BOUND) %>%
    suppressWarnings()
}

moldm <- raw_moldm("data/raw/db_20250616/novartis.csv")

# notes:
# S446I possible typo?
# haplotypes are falling out here but that's a super limited number
# mixeds where they are documented ...
# "C469C/Y", "C469F/Y", "M476M/I", "Y493Y/H", "R539R/T", ..........

moldm <- read.csv("data/clean/moldm_marcse_with_markers.csv")

message(paste("Number of studies:", length(unique(moldm$Title))))

# haplotypes? FIX review
message(paste("Number of rows indicating haplotypes:", 
              nrow(as.data.frame(moldm[grepl(",", moldm$Marker), c("Site.Name", "year", "Marker", "Present", "Tested")]))))

mutants <- moldm %>%
  filter(mutant) %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>%
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Present, Site.Name, Country, pubs) %>%
  ungroup() %>%
  suppressMessages()

moi <- moldm %>%
  filter(Marker == "A578S") %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>%
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Present, Site.Name, Country, pubs) %>%
  ungroup() %>%
  suppressMessages()
  

message(paste("Number of rows in mutant table:", nrow(mutants)))

# 578 takes off earlier
hist(mutants$year)
hist(moi$year)

plot(mutants$year, mutants$Present/mutants$Tested, xlab="Year", ylab="Prevalence")
plot(moi$year, moi$Present/moi$Tested, xlab="Year", ylab="Prevalence")

wildtypes <- moldm %>%
  filter(Marker == "k13 wildtype") %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), n=n(), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>% 
  # check how many simultaneous wildtype entries we have?
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country, Present, pubs) %>%
  ungroup() %>%
  suppressMessages()

# possible duplicates?
# FIX
wildtypes[which(wildtypes$Present/wildtypes$Tested > 1),
          c("Longitude", "Latitude",  "year", "Tested", "Site.Name",
            "Present", "pubs")]

wildtypes_to_add <- anti_join(
  # locations in `wildtypes` that do not occur in `mutants`,
  # paying attention to `Tested` but NOT to `pubs`
  wildtypes %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country),
  moi %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country)) %>%
  mutate(Present = 0) %>%
  suppressMessages()

nrow(wildtypes_to_add)

# write.csv(wildtypes_to_add,
#           "../k13_seafrica/data/clean/moldm_wildtypes_to_add.csv")

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
  labs(title = "Prevalence of 578S markers in Africa") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_bw() 
ggsave("../k13_ms/figs/A578S.png", height = 6, width = 9)
