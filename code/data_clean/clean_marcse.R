# tidy up MARCSE data in ../MARC_SEA_dashboard
# cloned from github: https://github.com/Stephanie-van-Wyk/MARC_SEA_dashboard.git
# (this mostly involves assigning PubMedIDs so that we can deduplicate during merge with moldm)

library(readxl)

# from MARC-SE GH README:
# markers <- "C580Y, P574L, R561H, P553L, I543T, R539T, Y493H, M476I, N458Y, 
#           F446I, C469Y, C469F, A675V, P441L, R622I, G625R, A578S, N537S" %>%
#   str_split(markers, ", ") %>% 
#   unlist()
# note 578 and 537 not included in WHO lists but included due to prevalence

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

marker_reference <- readxl::read_xlsx("compendium-of-molecular-markers-for-antimalarial-drug-resistance.xlsx",
                                      sheet = "Artemisinins (Pf)") %>%
  filter(grepl("Validated", Classification) | grepl("Candidate", Classification)) %>%
  rename(marker = `Alteration(s)`,
         status = Classification) %>%
  dplyr::select(marker, status)

marcse <- read_xlsx("Dashboard_k13_update_January_2026.xlsx") %>%
  mutate_at(c("Present", "Tested"), as.numeric) %>%
  dplyr::select(-c("...19": "...22")) %>%
  rename(Notes = `...18`,
         Site.Name.District.Country = `Site Name/District/Country`,
         Start.Year = `Start Year`,
         End.Year = `End Year`,
         Year.Published = `Year Published`,
         Prevalence... = `Prevalence (%)`,
         Site.Name = `Site Name`) %>%
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
  # not sure what happened here, but these titles and PMIDs aren't right ...
  mutate(Title = case_when(Title == "Comprehensive analysis of molecular markers linked to antimalarial drug resistance in Plasmodium falciparum in Northern. Northeastern and Eastern Uganda" & 
                             Authors == "Anthony Nuwa. Kevin Baker. Richard Kajubi. Chukwudi A Nnaji. Katherine Theiss-Nyland. Musa Odongo. Tonny Kyagulanyi. Jane Nabakooza. David Salandini. Victor Asua. Maureen Nakirunda. Christian Rassi. Damian Rutazaana. Richard Achuma. Patrick Sagaki. John Baptist Bwanika. Godfrey Magumba. Adoke Yeka. Sam Nsobya. Moses R Kamya. James Tibenderana. Jimmy Opigo"
                           ~ "Effectiveness of sulfadoxine-pyrimethamine plus amodiaquine and dihydroartemisinin-piperaquine for seasonal malaria chemoprevention in Uganda: a three-arm, open-label, non-inferiority and superiority, cluster-randomised, controlled trial",
                           Title == "Comprehensive analysis of molecular markers linked to antimalarial drug resistance in Plasmodium falciparum in Northern. Northeastern and Eastern Uganda" &
                             Authors == "Angwe MK. Mwebaza N. Nsobya SL. Vudriko P. Dralabu S. Omali D. Tumwebaze MA. Ocan M."
                           ~ "Day 3 parasitemia and Plasmodium falciparum Kelch 13 mutations among uncomplicated malaria patients treated with artemether-lumefantrine in Adjumani district, Uganda",
                           TRUE ~ Title)) %>%
  # and these authors weren't populated ...
  mutate(Authors = ifelse(is.na(Authors) & Title == "Comprehensive analysis of molecular markers linked to antimalarial drug resistance in Plasmodium falciparum in Northern. Northeastern and Eastern Uganda",
                          "Peter Olupot-Olupot  1   2 , George Paasi  3   4 , Thomas Katairo  5 , Jimmy Patrick Alunyo  3   4 , Alice Nakiyemba  6 , Gilbert Gilibrays Ocen  7 , Stephen Pande  8 , Florance Alaroker  9 , William Okiror  3   4 , Emmaluel Ocen  10 , Alex Oula  10 , Charles Benard Okalebo  3   4 , Ongodia Paul  3 , Denis Amorut  3 , Stephen Tukwasibwe  5 , Susan Nabadda Ndidde  11 , Isaac Sewanyana  11 , Samuel L Nsobya  5", Authors)) %>%
  # something has corrupted PubMedID field ........... or these are otherwise blank:
  mutate(PubMedID = case_when(
    Title == "A Novel Plasmodium falciparum Kelch13 A675T Mutation and High Levels of Chloroquine and Sulfadoxine-Pyrimethamine Resistance in Burundi" ~ "40666336",
    Title == "Antimalarial drug resistance and population structure of Plasmodium falciparum in Mozambique using genomic surveillance at health facilities (2021-2022)" ~ "40790052",
    Title == "Artemisinin Partial Resistance Mutations in Zanzibar and Tanzania Suggest Regional Spread and African Origins. 2023" ~ "40802860",
    Title == "Changes in susceptibility of Plasmodium falciparum to antimalarial drugs in Uganda over time: 2019-2024." ~ "40783405",
    Title == "Comprehensive analysis of molecular markers linked to antimalarial drug resistance in Plasmodium falciparum in Northern. Northeastern and Eastern Uganda" ~ "40514714",
    Title == "Effectiveness of sulfadoxine-pyrimethamine plus amodiaquine and dihydroartemisinin-piperaquine for seasonal malaria chemoprevention in Uganda: a three-arm, open-label, non-inferiority and superiority, cluster-randomised, controlled trial" ~ "39826559",
    Title == "Day 3 parasitemia and Plasmodium falciparum Kelch 13 mutations among uncomplicated malaria patients treated with artemether-lumefantrine in Adjumani district, Uganda" ~ "38837973",
    Title == "Detection of Twenty-Four Plasmodium Falciparum Kelch 13 Mutations Including C469Y. P553L. R561H. and A675V Across Kenya" ~ "10.2139/ssrn.5020665", # doi will have to do
    Title == "Detection of Low-Frequency Artemisinin Resistance Mutations C469Y. P553L and A675V in Asymptomatic Primary School Children in Kenya" ~ "39819451",
    Title == "Efficacies of artemether-lumefantrine. artesunate-amodiaquine. dihydroartemisinin-piperaquine. and artesunate-pyronaridine for the treatment of uncomplicated Plasmodium falciparum malaria in children aged 6 months to 10 years in Uganda: a randomised. open-label. phase 4 clinical trial." ~ "40845863",
    Title == "Efficacy and Safety of Artemether-Lumefantrine Against Uncomplicated Falciparum Malaria Infection in Tanzania. 2022: A Single-Arm Clinical Trial." ~ "39186698",
    Title == "Efficacy and Safety of Artesunate–Amodiaquine and Artemether–Lumefantrine for the Treatment of Uncomplicated Plasmodium falciparum Malaria in Madagascar. 2020" ~ "41187342",
    Title == "Efficacy of artesunate-amodiaquine and artemether-lumefantrine for uncomplicated Plasmodium falciparum malaria in Madagascar. 2022" ~ "36376921",
    Title == "Genomic Surveillance Reveals Clusters of Plasmodium falciparum Antimalarial Resistance Markers in Eswatini. a Low-Transmission Setting" ~ "10.1101/2025.07.30.25332463",
    Title == "Global assessment of partial artemisinin resistance: multicenter trial across Kenya. Peru. and Thailand in patients with uncomplicated Plasmodium falciparum malaria." ~ "40614930",
    Title == "High Prevalence of Molecular Markers Associated with Artemisinin. Sulphadoxine and Pyrimethamine Resistance in Northern Namibia" ~ "Already in moldm", #"40744004",
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
    Title == "Very low prevalence of validated kelch13 mutations and absence of hrp2/3 double gene deletions in South African malaria-eliminating districts (2022-2024)." ~ "Already in moldm", # "11998825",
    Title == "Molecular surveillance of Plasmodium falciparum drug resistance markers in clinical samples from Botswana." ~ "30350774",
    Title == "WHO Threats Map" ~ "WHO Threats Map",
    Title == "Unpublished. WHO Threats Map" ~ "WHO Threats Map",
    Title == "MIM Genomic surveillance of Plasmodium falciparum Drug Resistance Markers and Genetic Diversity in Mozambique" ~ "WHO Threats Map",
    Title == "Malaria update: Increase in frequency of Kelch 13 mutations found" ~ "Unpublished",
    Title == "https://www.nicd.ac.za/wp-content/uploads/2022/08/310822-NICD-Monthly-Communique-Aug-NW5.pdf" ~ "Unpublished",
    Title == "2025" ~ "Unpublished",
    Title == "Plasmodium falciparum Kelch-13 artemisinin partial resistance markers in Fort Portal, Western Uganda, 2024" ~ "40265952", # published
    TRUE ~ PubMedID)) %>%
  # present and tested are swapped for this study:
  mutate(Tested = ifelse(PubMedID == "39800719", Present, Tested),
         Present = ifelse(PubMedID == "39800719", 0, Present)) %>%
  suppressMessages()

message("Check over irregular PMIDs")
marcse %>%
  filter(str_length(PubMedID) != 8 & 
           PubMedID != "Unpublished" & 
           PubMedID != "WHO Threats Map") %>%
  dplyr::select(PubMedID, Title, Authors) %>%
  unique() %>%
  as.data.frame()

# marcse[,c("Marker", "Marker_Classification", "status")] %>% unique() %>% as.data.frame()
# A626S, A675V not picked up in LH's set
# R515K not picked up in SvW's set
# then there's also A578S
# there aren't any 537 classified as "Widespread other"


# there's one particularly weird pmid 9999971 but satisfied that it's not published ..


# there are a couple with weird PMIDs and Title == "WHO Threats Map"
# pmid_checks <- c(
#   "https://www.nicd.ac.za/wp-content/uploads/2022/08/310822-NICD-Monthly-Communique-Aug-NW5.pdf", "Unpublished",
#   "https://www.nicd.ac.za/wp-content/uploads/2021/08/NICD-Monthly-Communique%CC%81-August.pdf", "Unpublished",
#   "84.8", "39819451", # end up throwing these out - NA Tested
#   "96", "39819451",
#   "98.8", "39819451",
#   "98.3", "39819451",
#   "99.86", "Unpublished", # this one is 1 - prevalence?
#   "99.72", "Unpublished",
#   "Unpublished. WHO Threats Map", "WHO Threats Map",
#   "Unavailable", "Unpublished",
#   "92.33", "WHO Threats Map",
#   "94.1", "WHO Threats Map",
#   "95.3", "WHO Threats Map",
#   # these weird PMIDs are in the Surveyor set:
#   "3", "Unpublished",
#   "4", "Unpublished",
#   "2", "Unpublished") %>%
#   matrix(byrow = TRUE, ncol = 2) %>%
#   as.data.frame() %>%
#   setNames(c("from", "to"))

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

marcse %>%
  filter(PubMedID == "Unpublished") %>% 
  group_by(Title) %>% 
  summarise(n=n()) %>%
  as.data.frame()

# I've re-extracted study 40744006
reextractees <- bind_rows(
  read.csv("reextract/40744006_clean.csv"),
  read.csv("reextract/39802768_clean.csv")) %>%
  mutate(Marker = ifelse(Marker == "WT", "wildtype", Marker),
         year = round((Start.Year + End.Year) / 2, 0),
         across(c(PubMedID, Year.Published), ~ as.character(.x))) %>%
  dplyr::select(-c(Province)) 

marcse <- marcse %>%
  filter(PubMedID != "40744006" & PubMedID != "39802768") %>%
  bind_rows(reextractees) %>%
  left_join(marker_reference, by = join_by(Marker == marker))

message("check for sensibility: Testeds < Presents")
marcse %>%
  filter(Present > Tested) %>%
  dplyr::select(Present, Tested, Title, Authors, PubMedID, Marker)
  # group_by(Title, Authors, PubMedID) %>%
  # summarise(n = n()) %>%
  # as.data.frame()

# will have to remove the weird record from MTM,
tmp <- marcse %>%
  filter(!(Tested == 106 & Title == "WHO Threats Map" & Marker == "wildtype")) %>%
  filter()

# concerned about 39800719
marcse %>% filter(PubMedID == "39802768") %>% as.data.frame()






