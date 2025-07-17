# merge IDDO Surveyor data in data/raw/ with MARCSE data in ../MARC_SEA_dashboard
# cloned from github: https://github.com/Stephanie-van-Wyk/MARC_SEA_dashboard.git

# from MARC-SE GH README:
markers <- "C580Y, P574L, R561H, P553L, I543T, R539T, Y493H, M476I, N458Y, 
          F446I, C469Y, C469F, A675V, P441L, R622I, G625R, A578S, N537S" %>%
  str_split(markers, ", ") %>% 
  unlist()
# note 578 and 537 not included in WHO lists but included due to prevalence

marcse <- read_xlsx("../MARC_SEA_dashboard/k13_marcse_africa_GHR.xlsx") %>%
  suppressMessages() %>%
  rename(Notes = "...18",
         Site.Name.District.Country = `Site Name/District/Country`,
         Start.Year = `Start Year`,
         End.Year = `End Year`,
         Year.Published = `Year Published`,
         Prevalence... = `Prevalence (%)`) %>%
  mutate_at(c("Longitude", "Latitude", "Present", "Tested"), as.numeric) %>%
  suppressWarnings() %>%
# so the .csv is data that is exported to the dashboard - so all duplicates
  # bind_rows(read.csv("../MARC_SEA_dashboard/k13_marcse_africa.csv") %>%
  #             mutate_at(c("Longitude"), as.numeric))
  mutate(Marker = ifelse(Marker == "WT", "wildtype", Marker)) %>%
  left_join(marker_reference, by = join_by(Marker == marker))
head(marcse)
names(marcse)

# Bagamoyo Longitude: 38.900002? .. will need to reassign:

unique(marcse$PubMedID)
unique(marcse$Marker)
unique(marcse$Marker_Classification)
marcse[,c("Marker", "Marker_Classification", "status")] %>% unique() %>% as.data.frame()
# A626S, A675V not picked up in LH's set
# R515K not picked up in SvW's set
# then there's also A578S
# there aren't any 537 classified as "Widespread other"

# at this point, ya need to run the top of data_clean_k13.R
# ... I could put everything in the same place or put the functions somewhere on
# they're own but I'm going to keep things separate for now
moldm <- raw_moldm("data/raw/db_20250616/novartis.csv") %>%
  mutate(Marker = strip_marker)
head(moldm)

head(marker_reference)
