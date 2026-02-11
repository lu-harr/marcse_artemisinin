library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
library(terra)
library(iddoPal)
library(patchwork)
library(cowplot)

# marker_reference <- readxl::read_xlsx("data/marker_index.xlsx")
# from WHO markers compendium:
marker_reference <- readxl::read_xlsx("compendium-of-molecular-markers-for-antimalarial-drug-resistance.xlsx",
                                      sheet = "Artemisinins (Pf)") %>%
  filter(grepl("Validated", Classification) | grepl("Candidate", Classification)) %>%
  rename(marker = `Alteration(s)`,
         status = Classification) %>%
  dplyr::select(marker, status)


YEAR_LOWER_BOUND <- 2000
MIN_SAMPLE_SIZE <- 10

theme_set(theme_bw())

world <- ne_countries(scale="medium", returnclass = "sf")
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()


format_moldm_k13 <- function(path){
  # reads in fresh version of moldm, cleans up, and gives us some summary stats
  out <- read.csv(path) %>%
    mutate(across(c(Start.Year, End.Year, Present, Tested), as.numeric)) %>% 
    filter(!is.na(Start.Year) & !is.na(End.Year)) %>% # remove where both are NA
    filter(End.Year > 1960 & End.Year < 2500) %>%
    filter(Start.Year > 1960 & Start.Year < 2500) %>%
    # there definitely aren't any "Kelch 13" left in here are there? Nope
    filter(grepl("[k|K]13", Marker)) %>%
    filter(!(Longitude < -10 & Latitude < -10)) %>%
    filter(Continent == "Africa") %>%
    suppressWarnings()
  
  message(paste("Number of rows indicating double mutants:", 
                nrow(filter(out, grepl(",", Marker)))))
  
  # remove double mutants? - checked and they're not included as single mutants (see below)
  # there's probably a tidy way to do this but alas:
  double_mutants <- filter(out, grepl(",", Marker)) %>%
    separate_rows(Marker, sep = ",") %>%
    mutate(Marker = trimws(Marker))
  # I checked and there aren't any double mutants where both mutants are on the WHO list
  # So we don't end up counting anyone twice in aggregate model
  
  out <- filter(out, !grepl(",", Marker)) %>%
    # add disaggregated double mutants back in:
    bind_rows(double_mutants) %>%
    # strip "k13"
    mutate(strip_marker = gsub("[k|K]13 ", "", Marker)) %>% 
    # there's also a stray "\t" in there ...
    mutate(strip_marker = gsub("\\t", "", strip_marker)) %>%
    # handle these mixed infections carefully:
    # C469F/Y, C469Y/F, N537I/D, C469STOC/P
    # 469Y/F - count towards Y as I don't want to count it twice in aggregate model
    mutate(strip_marker = case_when(strip_marker == "C469F/Y" | strip_marker == "C469Y/F" ~ "C469Y",
                                    strip_marker == "N537I/D" ~ "N537I", # (this one's in the compendium)
                                    # "C469STOC/P" - this one will get filtered out anyway
                                    TRUE ~ strip_marker)) %>%
    # strip out mixed infections - count towards mutants
    mutate(strip_marker = gsub("./", "", strip_marker)) %>%
    # now link up with marker status
    left_join(marker_reference, 
              by = join_by(strip_marker == marker)) %>%
    mutate(year = round((Start.Year + End.Year) / 2, 0),
           # if one or the other is not complete, populate with the value we have:
           year = case_when(is.na(year) & !is.na(Start.Year) ~ Start.Year,
                            is.na(year) & !is.na(End.Year) ~ End.Year,
                            TRUE ~ year),
           Longitude = as.numeric(Longitude),
           Latitude = as.numeric(Latitude),
           mutant = !is.na(status)) %>%
    suppressWarnings()
  
  message(paste0("Number of studies before filtering early records: ", 
                 length(unique(out$Title))))
  message(paste0("Publication years: ", 
                 range(as.numeric(out$Year.Published), na.rm=TRUE) %>%
                   suppressWarnings()))
  message(paste("Earliest years:", min(out$year, na.rm=TRUE)))
  
  out <- filter(out, year >= YEAR_LOWER_BOUND)
  
  message(paste0("Number of studies after filtering early records: ", 
                 length(unique(out$Title))))
  # this is a little naive but it's the best estimate we're going to get:
  message(paste0("Number of people screened: ", 
                 out %>%
                   group_by(Longitude, Latitude, year, Title, Tested) %>% 
                   summarise(n = n()) %>%
                   ungroup() %>%
                   dplyr::select(Tested) %>%
                   sum() %>%
                   suppressMessages()))
  
  message(paste0("Or more conservatively: ", 
                 out %>%
                   group_by(Longitude, Latitude, year, PubMedID) %>% 
                   summarise(n = length(unique(Tested)), Tested = max(Tested)) %>%
                   ungroup() %>%
                   dplyr::select(Tested) %>%
                   sum() %>%
                   suppressMessages()))
  
  out
}