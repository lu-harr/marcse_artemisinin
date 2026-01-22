# merge IDDO Surveyor data in data/raw/ with MARCSE data in ../MARC_SEA_dashboard
# cloned from github: https://github.com/Stephanie-van-Wyk/MARC_SEA_dashboard.git

marker_reference <- readxl::read_xlsx("data/marker_index.xlsx")

# at this point, ya need to run the top of data_clean_k13.R
# ... I could put everything in the same place or put the functions somewhere on
# they're own but I'm going to keep things separate for now
moldm <- raw_moldm("data/raw/db_20260105/novartis.csv") %>%
  mutate(Marker = strip_marker)
# head(moldm)

# from MARC-SE GH README:
# markers <- "C580Y, P574L, R561H, P553L, I543T, R539T, Y493H, M476I, N458Y, 
#           F446I, C469Y, C469F, A675V, P441L, R622I, G625R, A578S, N537S" %>%
#   str_split(markers, ", ") %>% 
#   unlist()
# note 578 and 537 not included in WHO lists but included due to prevalence

library(readxl)

# marcse <- read_xlsx("../MARC_SEA_dashboard/k13_marcse_africa_GHR.xlsx") %>%
#   mutate_at(c("Present", "Tested"), as.numeric) %>%
#   bind_rows(read_xlsx("../MARC_SEA_dashboard/September_25_Lucy.xlsx") %>%
#               dplyr::select(-c("...19", "Site Name", "WIld-type (%)"))) %>%
#   rename(Notes = "...18",
#          Site.Name.District.Country = `Site Name/District/Country`,
#          Start.Year = `Start Year`,
#          End.Year = `End Year`,
#          Year.Published = `Year Published`,
#          Prevalence... = `Prevalence (%)`) %>%
#   mutate_at(c("Longitude", "Latitude"), as.numeric) %>%
#   mutate(year = round((Start.Year + End.Year) / 2, 0),
#          Year.Published = as.character(Year.Published) # clashing during bind_rows
#          ) %>%
#   suppressWarnings() %>%
# # so the .csv is data that is exported to the dashboard - so all duplicates
#   # bind_rows(read.csv("../MARC_SEA_dashboard/k13_marcse_africa.csv") %>%
#   #             mutate_at(c("Longitude"), as.numeric))
#   mutate(Marker = ifelse(Marker == "WT", "wildtype", Marker),
#          # in preparation for casting to numeric ..
#          Prevalence = gsub("%", "", Prevalence...), 
#          # in preparation for imputing Testeds
#          Prevalence = as.numeric(Prevalence),
#          # for one study, Tested was not provided but Present and Prev were
#          Tested = ifelse(Title == "Detection of Low-Frequency Artemisinin Resistance Mutations C469Y. P553L and A675V in Asymptomatic Primary School Children in Kenya",
#                          Present / (Prevalence / 100), Tested)) %>%
#   left_join(marker_reference, by = join_by(Marker == marker)) %>%
#   suppressMessages()

marcse <- read_xlsx("../MARC_SEA_dashboard/Dashboard_k13_update_January_2026.xlsx") %>%
  mutate_at(c("Present", "Tested"), as.numeric) %>%
  dplyr::select(-c("...19": "...22")) %>%
  rename(Notes = `...18`,
         Site.Name.District.Country = `Site Name/District/Country`,
         Start.Year = `Start Year`,
         End.Year = `End Year`,
         Year.Published = `Year Published`,
         Prevalence... = `Prevalence (%)`) %>%
  mutate_at(c("Longitude", "Latitude"), as.numeric) %>%
  mutate(year = round((Start.Year + End.Year) / 2, 0),
         Year.Published = as.character(Year.Published) # clashing during bind_rows
  ) %>%
  suppressWarnings() %>%
  # so the .csv is data that is exported to the dashboard - so all duplicates
  # bind_rows(read.csv("../MARC_SEA_dashboard/k13_marcse_africa.csv") %>%
  #             mutate_at(c("Longitude"), as.numeric))
  mutate(Marker = ifelse(Marker == "WT", "wildtype", Marker),
         # in preparation for casting to numeric ..
         Prevalence = gsub("%", "", Prevalence...), 
         # in preparation for imputing Testeds
         Prevalence = as.numeric(Prevalence),
         # for one study, Tested was not provided but Present and Prev were
         Tested = ifelse(Title == "Detection of Low-Frequency Artemisinin Resistance Mutations C469Y. P553L and A675V in Asymptomatic Primary School Children in Kenya",
                         Present / (Prevalence / 100), Tested),
         Present = ifelse(is.na(Present),
                          Tested * Prevalence,
                          Present)) %>%
  # something has corrupted PubMedID field ........... 
  mutate(PubMedID = case_when(
    Title == "A Novel Plasmodium falciparum Kelch13 A675T Mutation and High Levels of Chloroquine and Sulfadoxine-Pyrimethamine Resistance in Burundi" ~ "12262749",
    Title == "Antimalarial drug resistance and population structure of Plasmodium falciparum in Mozambique using genomic surveillance at health facilities (2021-2022)" ~ "12340125",
    Title == "Artemisinin Partial Resistance Mutations in Zanzibar and Tanzania Suggest Regional Spread and African Origins. 2023" ~ "40802860",
    Title == "Changes in susceptibility of Plasmodium falciparum to antimalarial drugs in Uganda over time: 2019-2024." ~ "12335539",
    Title == "Comprehensive analysis of molecular markers linked to antimalarial drug resistance in Plasmodium falciparum in Northern. Northeastern and Eastern Uganda" ~ "12164153",
    Title == "Detection of Twenty-Four Plasmodium Falciparum Kelch 13 Mutations Including C469Y. P553L. R561H. and A675V Across Kenya" ~ "10.2139/ssrn.5020665", # doi will have to do
    Title == "Efficacies of artemether-lumefantrine. artesunate-amodiaquine. dihydroartemisinin-piperaquine. and artesunate-pyronaridine for the treatment of uncomplicated Plasmodium falciparum malaria in children aged 6 months to 10 years in Uganda: a randomised. open-label. phase 4 clinical trial." ~ "40845863",
    Title == "Efficacy and Safety of Artemether-Lumefantrine Against Uncomplicated Falciparum Malaria Infection in Tanzania. 2022: A Single-Arm Clinical Trial." ~ "39186698",
    Title == "Efficacy and Safety of Artesunate–Amodiaquine and Artemether–Lumefantrine for the Treatment of Uncomplicated Plasmodium falciparum Malaria in Madagascar. 2020" ~ "41187342",
    Title == "Efficacy of artesunate-amodiaquine and artemether-lumefantrine for uncomplicated Plasmodium falciparum malaria in Madagascar. 2022" ~ "36376921",
    Title == "Genomic Surveillance Reveals Clusters of Plasmodium falciparum Antimalarial Resistance Markers in Eswatini. a Low-Transmission Setting" ~ "10.1101/2025.07.30.25332463",
    Title == "Global assessment of partial artemisinin resistance: multicenter trial across Kenya. Peru. and Thailand in patients with uncomplicated Plasmodium falciparum malaria." ~ "40614930",
    Title == "High Prevalence of Molecular Markers Associated with Artemisinin. Sulphadoxine and Pyrimethamine Resistance in Northern Namibia" ~ "40744004",
    Title == "KEMRI Report: TES report for Malawi July 2025" ~ "Unpublished",
    Title == "MIM Conference. Kawela M et al Preliminary Regional Results from GenE8" ~ "Unpublished",
    Title == "Malaria prevalence. transmission potential and efficacy of artemisinin-based combination therapy in the Kenyan Central highlands: a zone previously characterized as malaria-free." ~ "39800719",
    Title == "Pharmacometric evaluation of amodiaquine-sulfadoxine-pyrimethamine and dihydroartemisinin-piperaquine seasonal malaria chemoprevention in northern Uganda" ~ "41231725",
    Title == "Plasmodium falciparum Kelch-13 artemisinin partial resistance markers in Fort Portal. Western Uganda. 2024" ~ "40265952",
    Title == "Plasmodium falciparum genomic surveillance reveals a diversity of kelch13 mutations in Zambia" ~ "40744006",
    Title == "Prevalence of Plasmodium species in asymptomatic individuals in North-Eastern South Africa: 2018 - 2019." ~ "41378557",
    Title == "Prevalence of resistance markers of artemisinin. partner drugs. and sulfadoxine-pyrimethamine in Nanyumbu and Masasi Districts. Tanzania between 2020 and 2021." ~ "40938322",
    Title == "SAMEC SCAT Feb 2025" ~ "Unpublished",
    Title == "The E8-led Regional Malaria Molecular Surveillance Initiative: Successes. Challenges. and Opportunities" ~ "Unpublished",
    Title == "Very low prevalence of validated kelch13 mutations and absence of hrp2/3 double gene deletions in South African malaria-eliminating districts (2022-2024)." ~ "11998825",
    Title == "WHO Threats Map" ~ "WHO Threats Map",
    Title == "2025" ~ "Unpublished",
    TRUE ~ PubMedID)) %>%
  left_join(marker_reference, by = join_by(Marker == marker)) %>%
  suppressMessages()

# head(marcse)
# names(marcse)

# unique(marcse$PubMedID)
# unique(marcse$Marker)
# unique(marcse$Marker_Classification)

# marcse[,c("Marker", "Marker_Classification", "status")] %>% unique() %>% as.data.frame()
# A626S, A675V not picked up in LH's set
# R515K not picked up in SvW's set
# then there's also A578S
# there aren't any 537 classified as "Widespread other"

# let's check for weirdness here:
marcse %>%
  filter(!PubMedID %in% moldm$PubMedID) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()
marcse %>%
  filter(str_length(PubMedID) != 8) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()

# check for weirdness here:  
moldm %>%
  filter(!PubMedID %in% marcse$PubMedID) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()
moldm %>%
  filter(str_length(PubMedID) != 8) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()
# there's one particularly weird pmid 9999971 but satisfied that it's not published ..


# there are a couple with weird PMIDs and Title == "WHO Threats Map"
pmid_checks <- c(
  "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6283485/", "30350774",
  "https://www.nicd.ac.za/wp-content/uploads/2022/08/310822-NICD-Monthly-Communique-Aug-NW5.pdf", "Unpublished",
  "https://www.nicd.ac.za/wp-content/uploads/2021/08/NICD-Monthly-Communique%CC%81-August.pdf", "Unpublished",
  "84.8", "39819451", # end up throwing these out - NA Tested
  "96", "39819451",
  "98.8", "39819451",
  "98.3", "39819451",
  "99.86", "Unpublished", # this one is 1 - prevalence?
  "99.72", "Unpublished",
  "Unpublished. WHO Threats Map", "WHO Threats Map",
  "Unavailable", "Unpublished",
  "92.33", "WHO Threats Map",
  "94.1", "WHO Threats Map",
  "95.3", "WHO Threats Map",
  # these weird PMIDs are in the Surveyor set:
  "3", "Unpublished",
  "4", "Unpublished",
  "2", "Unpublished") %>%
  matrix(byrow = TRUE, ncol = 2) %>%
  as.data.frame() %>%
  setNames(c("from", "to"))

marcse %>%
  filter(str_length(PubMedID) != 8 & ! PubMedID %in% pmid_checks$from) %>%
  dplyr::select(Title, PubMedID) %>%
  unique() %>%
  as.data.frame()
                 
# filter(marcse, Title == "WHO Threats Map") %>%
#   dplyr::select(PubMedID) %>% as.data.frame()

# filter(moldm, PubMedID == "null") %>% dplyr::select("Title") %>% unique()

# "null":
# (Unpublished)                                                                                                                
# (Unpublished) Presence of k13 561H artemisinin resistance mutations in Plasmodium falciparum infections from Rwanda http://malariamatters.org/presence-of-k13-561h-artemisinin-resistance-mutations-in-plasmodium-falciparum-infections-from-rwanda/
# (38840941) Screening for antifolate and artemisinin resistance in Plasmodium falciparum clinical isolates from three hospitals of Eritrea
# (38440764) Increase of Plasmodium falciparum parasites carrying lumefantrine-tolerance molecular markers and lack of South East Asian pfk13 artemisinin-resistance mutations in samples collected from 2013 to 2016 in Côte d’Ivoire
# (Unpublished - thesis? Or 37349311) UNDERSTANDING RESIDUAL PLASMODIUM FALCIPARUM TRANSMISSION IN ZANZIBAR THROUGH MULTIPLEXED AMPLICON DEEP SEQUENCING
# (Preprint https://www.medrxiv.org/content/10.1101/2025.01.09.25320247v1.full-text) High Prevalence of Molecular Markers Associated with Artemisinin, Sulphadoxine and Pyrimethamine Resistance in Northern Namibia
# (38461239) Trends of Plasmodium falciparum molecular markers associated with resistance to artemisinins and reduced susceptibility to lumefantrine in Mainland Tanzania from 2016 to 2021
# (?) Antimalarial Drug Resistance Marker Prevalence Survey - 2016
# (doi: 10.61186/rabms.11.1.75) Molecular surveillance of artemisinin resistance-linked PFK13 gene polymorphisms in Adamawa State, Nigeria

# I've pasted this but forgotten what it belongs to "33864801"

# find "unpublished" rows
marcse %>%
  left_join(pmid_checks, by = join_by(PubMedID == from)) %>% 
  filter(PubMedID == "Unpublished") %>% 
  group_by(Title) %>% 
  summarise(n=n()) %>%
  as.data.frame()

marcse <- marcse %>%
  left_join(pmid_checks, by = join_by(PubMedID == from)) %>%
  # PubMedID == "0":
  mutate(PubMedID = case_when(Title == "SAMEC SCAT Feb 2025" ~ "Unpublished",
                              Title == "NMCP. Victor Asua" ~ "Unpublished",
                              Title == "https://www.nicd.ac.za/wp-content/uploads/2022/08/310822-NICD-Monthly-Communique-Aug-NW5.pdf" ~ "Unpublished", # difficult to know if this ended up elsewhere
                              Title == "WHO Threats Map" ~ "WHO Threats Map",
                              # Title == "Efficacy of artesunate-amodiaquine and artemether-lumefantrine for uncomplicated Plasmodium falciparum malaria in Madagascar. 2022" ~ "Unpublished", # see 34732201
                              # Title == "Plasmodium falciparum genomic surveillance reveals a diversity of kelch13 mutations in Zambia" ~ "40744006", # preprint: https://doi.org/10.1101/2025.02.19.25322554; AJTMH "ahead of print"?
                              # Title == "Detection of Twenty-Four Plasmodium Falciparum Kelch 13 Mutations Including C469Y. P553L. R561H. and A675V Across Kenya" ~ "Unpublished", # preprint: https://dx.doi.org/10.2139/ssrn.5020665
                              Title == "High Prevalence of Molecular Markers Associated with Artemisinin. Sulphadoxine and Pyrimethamine Resistance in Northern Namibia" ~ "Already in moldm", # preprint: https://doi.org/10.1101/2025.01.09.25320247
                              Title == "Antimalarial drug resistance and population structure of Plasmodium falciparum in Mozambique using genomic surveillance at health facilities (2021-2022)" ~ "40790052", # now published
                              Title == "Malaria update: Increase in frequency of Kelch 13 mutations found" ~ "Unpublished", # https://www.nicd.ac.za/wp-content/uploads/2022/08/Malaria-update.pdf
                              Title == "The E8-led Regional Malaria Molecular Surveillance Initiative: Successes. Challenges. and Opportunities" ~ "Unpublished", # appears to be a presentation given by Dr Jaishree Raman https://www.marcse-africa.org/news/blog/uniting-against-malaria-9th-southern-africa-research-conference
                              # Title == "Genomic Surveillance Reveals Clusters of Plasmodium falciparum Antimalarial Resistance Markers in Eswatini, a Low-Transmission Setting" ~ "Unpublished", # preprint - doi doesn't seem to be in moldm https://doi.org/10.1101/2025.07.30.25332463
                              # Title == "Comprehensive analysis of molecular markers linked to antimalarial drug resistance in Plasmodium falciparum in Northern, Northeastern and Eastern Uganda" ~ "40514714", # published
                              Title == "Plasmodium falciparum Kelch-13 artemisinin partial resistance markers in Fort Portal, Western Uganda, 2024" ~ "40265952", # published
                              # Title == "A Novel Plasmodium falciparum Kelch13 A675T Mutation and High Levels of Chloroquine and Sulfadoxine-Pyrimethamine Resistance in Burundi" ~ "40666336", # preprint
                              TRUE ~ PubMedID)) %>%
  mutate(PubMedID = case_when(!is.na(to) ~ to,
                              TRUE ~ PubMedID)) %>%
  # actually would like to tag these for later:
  mutate(PubMedID = case_when(PubMedID == "Unpublished" ~ "Unpublished_MARCSE",
                              TRUE ~ PubMedID),
         from = "marcse") %>%
 dplyr::select(-c(to))

#tmp %>% dplyr::select(PubMedID, to) %>% unique() %>% as.data.frame()

moldm <- moldm %>%
  left_join(pmid_checks, by = join_by(PubMedID == from)) %>%
  mutate(PubMedID = case_when(!is.na(to) ~ to,
                              TRUE ~ PubMedID)) %>%
  mutate(PubMedID = case_when(Title == "Antimalarial Drug Resistance Marker Prevalence Survey - 2016" ~ "Unpublished",
                              Title == "UNDERSTANDING RESIDUAL PLASMODIUM FALCIPARUM TRANSMISSION IN ZANZIBAR THROUGH MULTIPLEXED AMPLICON DEEP SEQUENCING" ~ "Unpublished",
                              Title == "Presence of k13 561H artemisinin resistance mutations in Plasmodium falciparum infections from Rwanda" ~ "Unpublished",
                              Title == "High Prevalence of Molecular Markers Associated with Artemisinin, Sulphadoxine and Pyrimethamine Resistance in Northern Namibia" ~ "Unpublished",
                              Title == "Screening for antifolate and artemisinin resistance in Plasmodium falciparum clinical isolates from three hospitals of Eritrea" ~ "38840941",
                              Title == "Increase of Plasmodium falciparum parasites carrying lumefantrine-tolerance molecular markers and lack of South East Asian pfk13 artemisinin-resistance mutations in samples collected from 2013 to 2016 in Côte d’Ivoire" ~ "38440764",
                              Title == "Trends of Plasmodium falciparum molecular markers associated with resistance to artemisinins and reduced susceptibility to lumefantrine in Mainland Tanzania from 2016 to 2021" ~ "38461239",
                              Title == "Investigation of Markers of Antimalarial Resistance During a Therapeutic Efficacy Study Conducted in Uganda, 2018–2019" ~ "Unpublished", # can't seem to find a trace of it .. maybe it was a conference presentation
                              TRUE ~ PubMedID),
         from = "moldm",
         Longitude = as.numeric(Longitude),
         Latitude = as.numeric(Latitude)) %>%
  dplyr::select(-c(to))

oldv <- moldm

to_add <- anti_join(marcse, moldm, by = join_by(PubMedID)) %>%
  filter(PubMedID != "Already in moldm")
message(paste0("Studies: ", length(unique(to_add$Title))))
message(paste0("Patients: ", to_add %>%
                 group_by(Longitude, Latitude, year, Title, Tested) %>%
                 summarise(n = n()) %>%
                 ungroup() %>%
                 dplyr::select(Tested) %>%
                 sum() %>%
                 suppressMessages()))
message(paste0("Or more conservatively: ", to_add %>%
                 group_by(Longitude, Latitude, year, Title) %>%
                 summarise(n = length(unique(Tested)), Tested = max(Tested)) %>%
                 ungroup() %>%
                 dplyr::select(Tested) %>%
                 sum() %>%
                 suppressMessages()))


moldm <- bind_rows(moldm,
                   anti_join(marcse, moldm, by = join_by(PubMedID))) %>%
  filter(PubMedID != "Already in moldm") %>% # a sneaky preprint snuck through
  mutate(mutant = !is.na(status) & status != "Not associated")

message("From here down it's pretty much as in the other cleaning script")

write.csv(moldm %>%
            filter(Present / Tested <= 1 & Tested > 5), 
          "data/clean/moldm_marcse_with_markers.csv")

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

message(paste("Number of rows in mutant table:", nrow(mutants)))
# added 1035 - 967 == 68 rows here

tmp <- moldm %>% 
  filter(mutant) %>% 
  group_by(Marker) %>%
  summarise(n_present = sum(Present), n_tested = sum(Tested)) %>%
  filter(n_present > 0) %>%
  #filter(npres >= 10) %>%
  full_join(marker_reference, join_by(Marker == marker)) %>%
  arrange(desc(n_present)) %>%
  filter(!is.na(n_present))
message("TF-associated mutations in dataset")
tmp %>% as.data.frame()



wildtypes <- moldm %>%
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

# poss duplicates?
# wildtypes[which(wildtypes$Present/wildtypes$Tested > 1),
#           c("Longitude", "Latitude",  "year", "Tested", "Site.Name",
#             "Present", "pubs")]

wildtypes_to_add <- anti_join(
  # locations in `wildtypes` that do not occur in `mutants`,
  # paying attention to `Tested` but NOT to `pubs`
  wildtypes %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country),
  mutants %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country)) %>%
  mutate(Present = 0) %>%
  suppressMessages()

write.csv(wildtypes_to_add,
          "data/clean/moldm_marcse_wildtypes_to_add.csv",
          row.names = FALSE)

message(paste("Number of rows of wildtypes to add:", nrow(wildtypes_to_add)))
# 536 - 413 == 123 added rows

with_wildtypes <- full_join(mutants, wildtypes_to_add) %>%
  filter(Tested > MIN_SAMPLE_SIZE) %>%
  suppressMessages()

# suspect this is present and tested around the wrong way? Check with Stephanie?
# (These are all from the same paper)
with_wildtypes %>%
  filter(Present / Tested > 1)
marcse %>%
  filter(Present / Tested > 1) %>%
  dplyr::select(PubMedID, Country, `Site Name`, Present, Tested, Marker, Prevalence, Title) %>%
  as.data.frame()

write.csv(with_wildtypes %>%
            dplyr::select(-c(pubs)) %>%
            filter(Present / Tested <= 1),
          "data/clean/moldm_marcse_k13_nomarker.csv",
          row.names = FALSE)

message(paste("Number of testees:", 
              with_wildtypes %>% 
                filter(Present / Tested <= 1) %>%
                group_by(year, Longitude, Latitude, Tested) %>% 
                summarise(n=n()) %>% 
                ungroup() %>%
                dplyr::select(Tested) %>% 
                sum()))
# 119666.02 ??? caused by imputing testeds from present/prevalence ??
# (Although could have sworn )
# used to be 82041

###############################################################################
with_wildtypes <- read.csv("data/clean/moldm_marcse_k13_nomarker.csv")
wildtypes_to_add <- read.csv("data/clean/moldm_marcse_wildtypes_to_add.csv")
moldm <- read.csv("data/clean/moldm_marcse_with_markers.csv")

plot(with_wildtypes$year, with_wildtypes$Present/with_wildtypes$Tested, 
     xlab="Year", ylab="Prevalence")

to_vis <- with_wildtypes %>% 
  mutate(year_bin = cut(year, 
                        breaks = c(min(year)-1, 2009, 2012, 2015, 2018, 2021, max(year)))) %>%
  # mutate(year_bin = case_when(year <= 2007 ~ "1995-2009",
  #                             year <= 2011 ~ "2010-2011",
  #                             year <= 2013 ~ "2012-2013",
  #                             year <= 2015 ~ "2014-2015",
  #                             year <= 2017 ~ "2016-2017",
  #                             year <= 2019 ~ "2018-2019",
  #                             year <= 2021 ~ "2020-2021",
  #                             TRUE ~ "2022-2024")) %>%
  arrange(Present/Tested, Tested) %>%
  drop_na(Longitude, Latitude)

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
  scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
  facet_wrap(~ year_bin, ncol=3) +
  labs(title = "Prevalence of k13 markers in Africa") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_grey() 
ggsave("figures/archive/abs_grey.png", height=3, width=5, scale=2)
# spooky warning here ...

markers <- moldm %>%
  filter(mutant) %>% # no filter on Tested
  group_by(year, Marker) %>%
  summarise(present = sum(Present)) %>%
  rename(marker = Marker)

# put test intensity in the background
bg <- moldm %>%
  group_by(Longitude, Latitude, year, PubMedID) %>%
  summarise(Tested = max(Tested)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(Tested = sum(Tested))

# I cannot state again how intensely easy this is in base graphics
markers_keep <- markers %>%
  filter(present > 5) %>%
  arrange(present) %>%
  ungroup() %>%
  dplyr::select(marker) %>%
  unique()

npal = 7
sty <- markers_keep %>%
  mutate(col = factor(rep(1:npal, 5)[1:length(marker)]),
         linet = factor(rep(c(1:2), each = npal)[1:length(marker)]))

df <- markers %>%
  filter(marker %in% markers_keep$marker) %>%
  left_join(sty, join_by(marker==marker)) %>%
  mutate(group = interaction(col, linet, sep = " / "))


# chop off <2010? fix colour palette some more?

# ggplot(markers %>%
#          filter(Marker %in% markers_keep$Marker)) +
#   geom_line(aes(x = year, y = present, colour = Marker)) +
#   labs(title="Most prevalent markers by year")

# Are there particular spatial trends between markers?
#markers_keep <- unique(markers$marker[markers$present > 25])
markers_keep <- markers %>% 
  ungroup() %>%
  group_by(marker) %>%
  summarise(n = sum(present)) %>%
  arrange(desc(n)) %>%
  dplyr::slice(1:5) %>%
  bind_rows(data.frame(marker = "Others", n=1))

wts <- wildtypes_to_add %>%
  dplyr::select(Longitude, Latitude, year, Tested) %>%
  filter(year > 2008) %>%
  merge(markers_keep, by = NULL) %>%
  #rename(y = marker) %>% # what was I doing with this again?
  rename(Marker = marker) %>%
  mutate(Present = 0)

markers_disagg <- moldm %>%
  mutate(Marker = gsub("[K|k]13 ", "", Marker)) %>%
  filter(mutant) %>%
  mutate(Marker = ifelse(Marker %in% markers_keep$marker,
                         Marker, "Others")) %>%
  filter(!(Marker == "Others" & Present == 0) & # don't want zeroes for all other markers, just WTs as BG
           Tested > 5) %>% # don't want prevalence == 1, tested < 5 points
  #filter(Marker %in% markers_keep$marker) %>%
  bind_rows(wts) %>% #  %>% mutate(Marker = "wt")
  mutate(year_bin = cut(year, breaks = c(min(year)-1, 2012, 2015, 2018, 2021, max(year)))) %>%
  filter(Longitude > -30) %>% # remove ocean point that's a bit far away
  mutate(Marker = factor(Marker, levels = markers_keep$marker))

iddoblue <- iddoPal::iddo_palettes$iddo[1]

# Combine color and linetype into aesthetics
# I would now like to control the order in which the lines are plotted
# It shouldn't be that hard
# I shouldn't have to spend longer than 5 mins on stackoverflow
# It doesn't matter that much
# This is how we end up with crap graphs !!!!!!
df <- df %>%
  filter(year >= 2005) %>%
  left_join(markers_keep, by = join_by(marker==marker)) %>%
  arrange(desc(n))

bg <- bg %>%
  filter(year >= 2005)

bg_scale <- 100
bg_col <- "grey65"

message("Caution: setting bg_scale and y limits manually")
p1 <- 
  ggplot() +
  geom_bar(data = bg, 
           aes(x = year, y = Tested / bg_scale), 
           stat = "identity", fill = "grey75") +
  geom_line(data = df, size = 0.7,
            aes(x = year, y = present, group = marker, color = marker, 
                linetype = marker)) +
  geom_point(data = df, 
             aes(x = year, y = present, group = marker, color = marker)) +
  labs(
    title = "(a) Detected mutations by year",
    color = "Marker",
    linetype = "Marker"
  ) +
  scale_color_manual(values = rep(c(iddoblue, "#c7047c", viridis(4), "#E37210"), 2)) +
  scale_linetype_manual(values = rep(1:2, each = 7)) +
  scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
                     limits = c(0, 400)) +
  theme_minimal() +
  xlab("Year") +
  ylab("Mutations detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(-0.3, "lines"),
        legend.box.margin = margin(50, 6, 6, 6),
        axis.text.y.right = element_text(color = bg_col),
        axis.title.y.right = element_text(color = bg_col),
        axis.ticks.y.right = element_line(color = bg_col)) +
  #legend.justification.right = "bottom") +
  scale_x_continuous(breaks = 2008:2025, limits = c(2008, 2025)) 
p1

p2 <- ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(markers_disagg, Present == 0),
             mapping = aes(x = Longitude, y = Latitude,
                           size = Tested, col = "grey50"),
             fill = "grey70",  pch = 21, alpha = 0.5, stroke = 0.2) +
  geom_point(data = markers_disagg %>% 
               filter(Present > 0) %>%
               arrange(Present / Tested), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2, alpha = 0.5) +
  scale_color_manual(name = "", values = c("grey30"), labels=c("Absence"), guide = "none") +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt") +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt") +
  # allows labelling of rows and columns:
  facet_grid(rows=vars(year_bin), cols=vars(Marker)) +
  #facet_wrap(~ year_bin + Marker, drop=FALSE) +
  labs(title = "(b) Prevalence of k13 markers in Africa") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_bw() +
  theme(plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(10, "lines"),
        panel.spacing = unit(0, "lines")) +
  guides(fill = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top"))

p2

# p2 +
#    geom_rect(data = data.frame(xmin = 60, xmax = 80, ymin = 0, ymax = 40, year_bin = "(2015,2018]",
#                                Marker = "R622I"), 
#              aes(xmin=xmin, xmax=xmax, ymin=ymin,  ymax=ymax), 
#              colour="grey10", fill="grey85", linewidth=0.3)

library(grid)
library(cowplot)
# rect <- rectGrob(
#   x = 0.922,
#   y = 0.27,
#   width = unit(0.2, "in"),
#   height = unit(0.2, "in"),
#   gp = gpar(fill = "grey70", colour="grey50", alpha = 0.5)
# )

rect <- rectGrob(
  # x = 0.915,
  # y = 0.27,
  x = 0.4, 0.05,
  width = unit(0.2, "in"),
  height = unit(0.2, "in"),
  gp = gpar(fill = "grey70", colour="grey50", alpha = 0.5)
)

p2 <- ggdraw(p2) +
  draw_grob(rect) + 
  # draw_label("Absence", x = 0.96, y = 0.27, size = 10)
  draw_label("Absence", x = 0.45, y = 0.05, size = 10)
#theme(#strip.background = element_blank(),
#strip.text.x = element_blank(),
#strip.placement = "outside") +
#ggsave("figs/markers_disagg.png", height = 4, width = 5, scale = 2)

plot_grid(p1, p2, ncol = 1, rel_heights = c(0.6, 2))

ggsave("figures/markers_disagg_marcse.png", height = 6, width = 5, scale = 2)



################################################################################

ggplot() +
  geom_bar(data = bg, aes(x = year, y = Tested / bg_scale), 
           stat = "identity", fill = "grey75") +
  geom_line(data = df, size = 0.7,
            aes(x = year, y = present, group = marker, color = marker, linetype = marker)) +
  geom_point(data = df, 
             aes(x = year, y = present, group = marker, color = marker)) +
  labs(
    title = "Detected mutations by year - IDDO Surveyor + MARCSE-Africa data",
    color = "Marker",
    linetype = "Marker"
  ) +
  scale_color_manual(values = rep(c(viridis(4), "#E37210", iddoblue, "#c7047c"), 2)) +
  scale_linetype_manual(values = rep(1:2, each = 7)) +
  scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
                     limits = c(0, 325)) +
  theme_minimal() +
  xlab("Year") +
  ylab("Mutations detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(-0.3, "lines"),
        legend.box.margin = margin(50, 6, 6, 6),
        axis.text.y.right = element_text(color = bg_col),
        axis.title.y.right = element_text(color = bg_col),
        axis.ticks.y.right = element_line(color = bg_col)) +
  #legend.justification.right = "bottom") +
  scale_x_continuous(breaks = seq(2005, 2025, 2))
ggsave("~/Desktop/presentations/MARCSE/moldm_marcse_k13_time.png", scale = 0.8, width = 10, height = 5)



dat <- read.csv("../k13_seafrica/data/clean/moldm_marcse_k13_nomarker.csv") %>%
  filter(year >= 2000) %>%
  mutate(year_bin = cut(year, breaks = seq(2000, 2025, 5))) %>%
  arrange(Present/Tested)

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(dat, Present == 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, col = "grey50"),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  geom_point(data = filter(dat, Present > 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_color_manual(name = "Absence", values = c("grey30"), labels=c("")) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt") +
  scale_size_continuous(name = "Tested", range = c(0.5, 6), trans = "sqrt") +
  #facet_wrap(~ year_bin, ncol=3) +
  labs(title = "Prevalence of Kelch 13 markers - with unpublished records") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_grey() 
ggsave("~/Desktop/presentations/MARCSE/moldm_marcse_k13.png", height = 6, width = 6)


################################################################################
# 
# afr_with_centroids <- afr %>%
#   st_centroid() %>%
#   st_coordinates() %>%
#   bind_cols(afr) %>%
#   mutate(region = case_when(Y < -10 ~ "south",
#                             X > 22 & Y > -10 & Y < 22 ~ "east",
#                             X < 22 & Y > -10 & Y < 22 ~ "west")) %>%
#   st_as_sf()
# 
# safr <- afr_with_centroids %>%
#   filter(Y < -10) %>%
#   st_as_sf()
# 
# eafr <- afr_with_centroids %>%
#   filter(X > 22 & Y > -10 & Y < 22) %>%
#   st_as_sf()
# 
# wafr <- afr_with_centroids %>%
#   filter(X < 22 & Y > -10 & Y < 22) %>%
#   st_as_sf()
# 
# dat_tagged <- markers_disagg %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(afr))
# 
# dat_tagged <- dat_tagged %>%
#   mutate(eafr = st_intersects(., eafr) %>% lengths(),
#          wafr = st_intersects(., wafr) %>% lengths(),
#          safr = st_intersects(., safr) %>% lengths(),
#          region = case_when(eafr == 1 ~ "east",
#                             wafr == 1 ~ "west",
#                             safr == 1 ~ "south",
#                             TRUE ~ NA)) %>%
#   mutate(year_bin = cut(year, breaks = c(min(year) - 1, 2015, 2018, 2021, 2024)))
# 
# p <- ggplot() + 
#   geom_sf(data = afr_with_centroids, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0),
#           mapping = aes(size = Tested, col = "grey50"),
#           fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#             filter(Present > 0) %>%
#             arrange(Present / Tested), 
#           mapping = aes(size = Tested, fill = Present / Tested),
#           col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt",
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                   max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   # allows labelling of rows and columns:
#   facet_grid(region ~ year_bin, scales = "free") +
#   #labs(title = "(b) Prevalence of k13 markers in Africa") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# 
# p
# 
# p_safr <- ggplot() + 
#   geom_sf(data = safr, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0 & safr == 1),
#              mapping = aes(size = Tested, col = "grey50"),
#              fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#                filter(Present > 0 & safr == 1) %>%
#                arrange(Present / Tested), 
#              mapping = aes(size = Tested, fill = Present / Tested),
#              col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt",
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                   max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   # allows labelling of rows and columns:
#   #facet_grid(rows=vars(year_bin), cols=vars(Marker)) +
#   #facet_wrap(~ year_bin + Marker, drop=FALSE) +
#   facet_wrap(~year_bin, ncol = 1) +
#   #labs(title = "(b) Prevalence of k13 markers in Africa") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# 
# p_wafr <- ggplot() + 
#   geom_sf(data = wafr, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0 & wafr == 1),
#           mapping = aes(size = Tested, col = "grey50"),
#           fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#             filter(Present > 0 & wafr == 1) %>%
#             arrange(Present / Tested), 
#           mapping = aes(size = Tested, fill = Present / Tested),
#           col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt", 
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                 max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   facet_wrap(~year_bin, ncol = 1) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# 
# p_eafr <- ggplot() + 
#   geom_sf(data = eafr, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0 & eafr == 1),
#           mapping = aes(size = Tested, col = "grey50"),
#           fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#             filter(Present > 0 & eafr == 1) %>%
#             arrange(Present / Tested), 
#           mapping = aes(size = Tested, fill = Present / Tested),
#           col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt", 
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                   max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   facet_wrap(~year_bin, ncol = 1) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# p_eafr
# 
# plot_grid(p_wafr, p_eafr, p_safr, nrow = 1)
