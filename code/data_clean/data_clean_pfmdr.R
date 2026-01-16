library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)

YEAR_LOWER_BOUND <- 2000
MIN_SAMPLE_SIZE <- 10

world <- ne_countries(scale="medium", returnclass = "sf")
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

moldm <- #read.csv("data/raw/db_20250616/novartis.csv") 
  read.csv("data/raw/db_20260105/novartis.csv") %>%
  mutate(across(c(Start.Year, End.Year, Present, Tested, Longitude, Latitude), 
                as.numeric)) %>% 
  filter(!is.na(Start.Year) & !is.na(End.Year)) %>% # remove where both are NA
  filter(End.Year > 1960 & End.Year < 2500) %>%
  filter(Start.Year > 1960 & Start.Year < 2500) %>%
  mutate(year = round((Start.Year + End.Year) / 2, 0),
         # if one or the other is not complete, populate with the value we have:
         year = case_when(is.na(year) & !is.na(Start.Year) ~ Start.Year,
                          is.na(year) & !is.na(End.Year) ~ End.Year,
                          TRUE ~ year)) %>%
  filter(!(Longitude < 1 & Longitude > 0 & Latitude > 0 & Latitude < 1)) %>%
  # this one is in Tanzania and clearly the Longitude is just lost:
  filter(!(Longitude < -10 & Longitude > -11 & Latitude > -11 & Latitude < -10)) %>%
  filter(Continent == "Africa") %>%
  suppressWarnings()

message("filtering some weird coords")

pfmdr <- moldm %>%
  filter(grepl("pfmdr1", Marker)) %>%
  mutate(Marker = gsub(pattern = "pfmdr1 ", replacement = "", Marker)) %>%
  dplyr::select(-c(Continent, Site.Name.District.Country, Prevalence....,
                   Prevalence.....incl.mixed, Contributer, Journal, Estimated.Year,
                   Estimated.Prev, Estimated.Location))
# drop where Display.On.Surveyor is false? Ask Sabina about how this works for pfmdr1



str_to_hap <- data.frame(stri = unique(pfmdr$Marker)) %>%
  mutate(tri = case_when(str_length(stri) == 5 & 
                           grepl("^[NY][FY]..[DY]$", stri) ~ 
                           paste0(substr(stri, 1, 2), substr(stri, 5, 5)),
                         TRUE ~ stri),
         `86` = case_when(grepl("86Y", tri) ~ 1,
                          grepl("86N/Y", tri) ~ 1,
                          grepl("N86", tri) ~ 0,
                          # where 86-184-1246 reported:
                          str_length(tri) == 3 & substr(tri, 1, 1) == "N" ~ 0,
                          str_length(tri) == 3 & substr(tri, 1, 1) == "Y" ~ 1,
                          # where 86-184-1034-1042-1246 reported:
                          # str_length(tri) == 5 & substr(tri, 1, 1) == "N" & 
                          #   tri != "N1042" ~ 0,
                          # str_length(tri) == 5 & substr(tri, 1, 1) == "Y" ~ 1,
                          TRUE ~ NA),
         `184` = case_when(grepl("184F", tri) ~ 1,
                           grepl("184Y/F", tri) ~ 1,
                           grepl("Y184", tri) ~ 0,
                           # where 86-184-1246 reported:
                           str_length(tri) == 3 & substr(tri, 2, 2) == "F" ~ 1,
                           str_length(tri) == 3 & substr(tri, 2, 2) == "Y" ~ 0,
                           # where 86-184-1034-1042-1246 reported:
                           # str_length(tri) == 5 & substr(tri, 2, 2) == "F" ~ 1,
                           # str_length(tri) == 5 & substr(tri, 2, 2) == "Y" ~ 0,
                           TRUE ~ NA),
         `1246` = case_when(grepl("1246Y", tri) ~ 1,
                            grepl("1246D/Y", tri) ~ 1,
                            grepl("D1246", tri) | grepl("1246D", tri) ~ 0,
                            # where 86-184-1246 reported:
                            str_length(tri) == 3 & substr(tri, 3, 3) == "D" ~ 0,
                            str_length(tri) == 3 & substr(tri, 3, 3) == "Y" & 
                              tri != "86Y" ~ 1,
                            # where 86-184-1034-1042-1246 reported:
                            # str_length(tri) == 5 & substr(tri, 5, 5) == "D" & 
                            #   tri != "1042D" ~ 0,
                            # str_length(tri) == 5 & substr(tri, 5, 5) == "Y" & 
                            #   !tri %in% c("1226Y", "1034Y", "86N/Y") ~ 1,
                            TRUE ~ NA))

str_to_hap %>%
  group_by(tri, `86`, `184`, `1246`) %>%
  summarise(n = n()) %>%
  as.data.frame()

str_to_hap %>%
  filter(is.na(`86`) & is.na(`184`) & is.na(`1246`))
# what to do with FD, FYD, 86?

pfmdr <- left_join(pfmdr, str_to_hap, by = join_by(Marker == stri)) %>%
  filter(!is.na(`86`) | !is.na(`184`) | !is.na(`1246`)) %>%
  mutate(n_loci = 3 - rowSums(is.na(dplyr::select(., c(`86`, `184`, `1246`)))))

nrow(pfmdr)

# these are rows where all three loci sequenced:
pfmdr %>%
  filter(!is.na(`86`) & !is.na(`184`) & !is.na(`1246`)) %>%
  nrow()

tmp <- pfmdr %>%
  filter(!is.na(`86`) & !is.na(`184`) & !is.na(`1246`)) %>%
  group_by(Longitude, Latitude, year, PubMedID, Tested) %>%
  summarise(n = n())
sum(tmp$Tested)

tmp <- pfmdr %>%
  filter(!(!is.na(`86`) & !is.na(`184`) & !is.na(`1246`))) %>%
  group_by(Longitude, Latitude, year, PubMedID, Tested) %>%
  summarise(n = n())
sum(tmp$Tested)

message("Here's the numbers for Results")
tmp <- pfmdr %>%
  filter(!is.na(`86`) | !is.na(`184`) | !is.na(`1246`)) %>%
  group_by(Longitude, Latitude, year, Title, Year.Published) %>%
  summarise(n = n(), Tested = max(Tested))
sum(tmp$Tested)

length(unique(tmp$Title))
length(setdiff(unique(tmp$Title), unique(crt$Title)))
range(tmp$year)
range(as.numeric(tmp$Year.Published), na.rm=TRUE)

# run assumptions past Sabina
# group into studies and review haplotypes/evidence of surveillance at all three loci
# I'm prepared to drop records if there's enough data with clear evidence of 
# surveillance at all three?
# how many reports are there of haplotypes?

par(mfrow = c(3,1), oma=c(0,0,2,0), mar = rep(4, 4))

pfmdr %>%
  filter(!is.na(`86`)) %>%
  # make this Marker to see the whole thing:
  dplyr::select(tri) %>%
  unlist() %>%
  factor() %>%
  table() %>%
  barplot(las=2, main="N86Y")

pfmdr %>%
  filter(!is.na(`184`)) %>%
  dplyr::select(tri) %>%
  unlist() %>%
  factor() %>%
  table() %>%
  barplot(las=2, main="Y184F")

pfmdr %>%
  filter(!is.na(`1246`)) %>%
  dplyr::select(tri) %>%
  unlist() %>%
  factor() %>%
  table() %>%
  barplot(las=2, main="D1246Y")
mtext("Reported haplotype frequencies", outer = TRUE)

# latent haplotypes:
latent <- expand.grid(c("Y", "N"), c("F", "Y"), c("Y", "D")) %>%
  mutate(latent = paste0(Var1, Var2, Var3))

# realised haplotypes (in the dataset):
realised <- pfmdr %>%
  group_by(`86`, `184`, `1246`) %>%
  summarise(n = n()) %>%
  mutate(a = ifelse(is.na(`86`), "X", ifelse(`86`, "Y", "N")),
         b = ifelse(is.na(`184`), "X", ifelse(`184`, "F", "Y")),
         c = ifelse(is.na(`1246`), "X", ifelse(`1246`, "Y", "D")),
         realised = paste0(a, b, c)) %>%
  dplyr::select(-c(a,b,c))

find_matches <- function(realised, cands, blank = "X"){
  # of candidate latent haplotypes, which match realised?
  for (loc in 1:str_length(realised)){
    if (substr(realised, loc, loc) != blank){
      cands <- cands[-which(substr(cands, loc, loc) != substr(realised, loc, loc))]
    }
  }
  cands
}

config <- sapply(realised$realised, function(x){
  out <- find_matches(x, latent$latent)
  ifelse(latent$latent %in% out, 1, 0)
}) %>%
  as.data.frame() %>%
  set_rownames(latent$latent)

pfmdr <- pfmdr %>%
  inner_join(realised, by = join_by(`86`, `184`, `1246`))

sing_trips <- pfmdr %>%
  group_by(Longitude, Latitude, year, Tested, PubMedID, n_loci) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Longitude, Latitude, year, n_loci) %>%
  summarise(Tested = sum(Tested), n = n()) %>%
  filter(n_loci != 2) # there's only one of these .. K at 86 didn't get counted

ggplot() + 
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_point(data = sing_trips %>%
               mutate(n_loci = ifelse(n_loci == 1, "Single locus", "Triple loci")) %>%
               arrange(desc(Tested)), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, fill = year),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_fill_viridis_c(name = "Year") +
  scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
  facet_wrap(~ n_loci) +
  #labs(title = "Sequencing at single and triple locus") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()
ggsave("figures/mdr_single_triple.png", height=4, width=8)

single_opps <- list("NXX"="YXX", "XXD"="XXY", "XYX"="XFX")
tmp <- names(single_opps)
names(tmp) <- single_opps

backwards <- names(single_opps)  # (wildtypes)
forwards <- names(tmp)  # (mutants) 
single_opps <- c(single_opps, tmp)

pfmdr <- pfmdr %>%
  filter(Tested >= MIN_SAMPLE_SIZE) %>%
  group_by(Longitude, Latitude, year, PubMedID, Tested, realised, `86`, `184`, `1246`, n_loci) %>%
  summarise(pres = paste(Present, collapse = ","), 
            n_tri = length(unique(tri)),
            n = n(),
            Present = ifelse(n > 1 & n_tri == 1, 
                             ifelse(length(unique(Mixed.Included)) > 1, # fingers crossed here?
                                    min(Present),
                                    max(Present)), 
                             sum(Present)),
            tri = paste(tri, collapse = ","), 
            marker = paste(Marker, collapse = ",")) %>%
  dplyr::select(-c(n, n_tri)) %>%
  ungroup()

# true duplicate:
pfmdr %>% filter(PubMedID == "26597254")
# same data entered twice? Not assigned different uniq_publication_ids?
# site: Cape Coast, 2008, 86Y
pfmdr %>% filter(PubMedID == "21288822")
# quintuple sequencing - not duplicates. can't currently get at this one
pfmdr %>% filter(PubMedID == "11751128")
# genuine mistake? nope - it's a year rounding problem maybe?
pfmdr %>% filter(PubMedID == "30559133" & year == "2016" & Tested == "43")

# otherwise I think I'm happy with this ...

# now work on opposites ...
single_loc <- pfmdr %>%
  filter(n_loci == 1) %>%
  mutate(`86` = !is.na(`86`),
         `184` = !is.na(`184`), 
         `1246` = !is.na(`1246`)) %>%
  group_by(Longitude, Latitude, year, PubMedID, Tested, `86`, `184`,`1246`) %>%
  summarise(n = n(),
            Present = ifelse(n == 2,
                             ifelse(realised[1] %in% forwards, Present[1], Present[2]),
                             ifelse(realised %in% forwards, Present, Tested - Present)),
            realised = ifelse(n == 2,
                              ifelse(realised[1] %in% forwards, realised[1], realised[2]),
                              ifelse(realised %in% forwards, 
                                     realised, unlist(single_opps[realised]))))

# given haplotype (string), number present (int) and configuration mat (rownames
# are full haplotypes; colnames are partial haplotypes), disaggregate "full"
# triple haplotype records into contributing partial haplotypes
disagg_trip <- function(haplot, present, config_mat){
  out <- config[1,] * 0
  for (hap in 1:length(haplot)){
    out <- out + config[haplot[hap],] * present[hap]
  }
  out
}

trips <- pfmdr %>%
  filter(n_loci == 3) %>%
  group_by(Longitude, Latitude, year, PubMedID, Tested) %>%
  summarise(disagg_trip(realised, Present, config)) %>%
  # keep mutant markers only:
  dplyr::select(Longitude, Latitude, year, PubMedID, Tested, YXX, XFX, XXY) %>%
  pivot_longer(c(YXX, XFX, XXY), names_to = "realised", values_to = "Present")

# do anti_join on lat, lon, year, realised, pub, (not tested)
to_add_from_trips <- anti_join(trips,
                               single_loc %>%
                                 ungroup() %>%
                                 dplyr::select(Longitude, Latitude, year, 
                                               realised, PubMedID)) %>%
  left_join(single_loc %>% 
              ungroup() %>% 
              dplyr::select(realised, `86`, `184`, `1246`) %>% 
              unique())

# huzzah! squeezed out a further 90 records out
single_loc <- single_loc %>% 
  bind_rows(to_add_from_trips) %>%
  mutate(loc = case_when(`86` ~ "N86Y",
                         `184` ~ "Y184F",
                         `1246` ~ "D1246Y",
                         TRUE ~ NA),
         year_bin = cut(year, breaks = c(1974, 2008, 2013, 2017, 2020, 2024))) %>%
  filter(Present/Tested <= 1 & Present/Tested >= 0) %>%
  arrange(desc(Tested))

write.csv(single_loc, "data/clean/pfmdr_single_locus.csv", row.names = FALSE)

write.csv(single_loc %>% filter(`86`), "data/clean/pfmdr_single_mdr86.csv", row.names = FALSE)
write.csv(single_loc %>% filter(`184`), "data/clean/pfmdr_single_mdr184.csv", row.names = FALSE)
write.csv(single_loc %>% filter(`1246`), "data/clean/pfmdr_single_mdr1246.csv", row.names = FALSE)

library(iddoPal)

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = single_loc %>% arrange(desc(Tested)), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, fill = Present/Tested),
             col = "grey50", pch=21, stroke = 0.2, alpha = 0.5) +
  scale_fill_gradientn(colors = iddoPal::iddo_palettes$BlGyRd, 
                       "",
                       breaks = c(0, 0.5, 1), 
                       labels = c("0  (all WT)", "0.5", "1  (all mutant)"),
                       limits = c(0,1)) +
  scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
  facet_grid(year_bin ~ loc) +
  labs(title = "Pfmdr1 84-186-1246: Surveillance at single loci") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_grey()
# ggsave("figures/pfmdr_single_locus.png", scale = 1.7, height = 6, width = 5)

plot(single_loc$year, single_loc$Present/single_loc$Tested, col = factor(single_loc$loc))

# problem here:
which(single_loc$Present / single_loc$Tested > 1)

tmp <- pfmdr %>%
  filter(n_loci == 3) %>%
  dplyr::select(Longitude, Latitude, year, PubMedID, Tested, Present, realised, 
                `86`, `184`, `1246`) %>%
  bind_rows(single_loc %>%
              dplyr::select(-c(n, loc, year_bin)))

barplot(table(tmp$realised))

tmp2 <- tmp %>%
  group_by(realised) %>%
  summarise(n = sum(Tested))

barplot(names.arg = tmp2$realised, height = tmp2$n)

write.csv(tmp, "data/clean/pfmdr_trip.csv", row.names = FALSE)

###############################################################################
# combine with crt to have a "markers for partner drug resistance" figure
# (need to have already run data_clean_crt - there's no reason for them to be
# separate other than that I wasn't thinking about mdr when I started thinking
# about crt)
crt <- read.csv("data/clean/moldm_crt76.csv") %>%
  mutate(loc = "Pfcrt K76T",
         year_bin = cut(year, breaks = c(1974, 2008, 2013, 2017, 2020, 2024)))

single_loc <- read.csv("data/clean/pfmdr_single_locus.csv")

message("Summary stats for results")

partners <- bind_rows(crt,
                      single_loc %>%
                        mutate(loc = paste("Pfmdr1", loc))) %>%
  mutate(loc = factor(loc, levels = c("Pfcrt K76T", "Pfmdr1 N86Y", "Pfmdr1 Y184F", "Pfmdr1 D1246Y"))) %>%
  filter(year > YEAR_LOWER_BOUND) %>%
  mutate(year_bin = cut(year, breaks = c(min(year) - 1, 2008, 2012, 2016, 2020, max(year))))

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = partners %>% arrange(desc(Tested)), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, fill = Present/Tested),
             col = "grey50", pch=21, stroke = 0.2, alpha = 0.5) +
  scale_fill_gradientn(colors = iddoPal::iddo_palettes$BlGyRd, 
                       "Prevalence",
                       breaks = c(0, 0.5, 1), 
                       labels = c("0\n(all wildtype)", "0.5", "1\n(all mutant)"),
                       limits = c(0,1)) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt",
                        breaks = c(10, 100, 1000, 3000)) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  facet_grid(year_bin ~ loc) +
  # labs(#title = "Pfmdr1 84-186-1246: Surveillance at single loci",
  #      xlab = "Longitude", ylab = "Latitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(4, "lines"),
        panel.spacing = unit(0, "lines")) +
  guides(fill = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top"))
#theme(title = element_blank())
ggsave("figures/crt_pfmdr_data.png", height = 9, width = 7)


partners <- partners %>%
  mutate(year_bin = cut(year, breaks = c(min(year) - 1, seq(2005, 2025, length.out = 5))))

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = partners %>%
               arrange(desc(Tested)), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, fill = Present/Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_fill_gradientn(colors = iddoPal::iddo_palettes$BlGyRd, 
                       "",
                       breaks = c(0, 0.5, 1), 
                       labels = c("0  (all wildtype)", "0.5", "1  (all mutant)"),
                       limits = c(0,1)) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt") +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  facet_grid(loc ~ year_bin, switch = "y") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(title = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.spacing.x = unit(0, "mm"),
        panel.spacing.y = unit(2, "mm"))
ggsave("figures/crt_pfmdr_data_short.png", scale = 1.7, height = 3, width = 4)

