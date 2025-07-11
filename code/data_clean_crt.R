# this is in fact my script for crt76
library(dplyr)
library(tidyr)
library(iddoPal)
# introduced bug to summarisation step: conflicting presences/apparent duplicates

# YEAR_LOWER_BOUND <- 2000
MIN_SAMPLE_SIZE <- 10

# here is my DIY set of verified/associated markers

# pfcrt 76: need 76T, K76, 76K/T, 
# ... not sure what to do about haplotypes? what about 74 and 75

crt <- read.csv("data/raw/db_20250616/novartis.csv") %>%
  mutate(Start.Year = as.numeric(Start.Year),
         End.Year = as.numeric(End.Year),
         Present = as.numeric(Present),
         Tested = as.numeric(Tested),
         Longitude = as.numeric(Longitude),
         Latitude = as.numeric(Latitude)) %>%
  filter(!is.na(Start.Year) & !is.na(End.Year)) %>% # remove where both are NA
  filter(End.Year > 1960 & End.Year < 2500) %>%
  filter(Start.Year > 1960 & Start.Year < 2500) %>%
  #filter(grepl("crt", Marker)) %>%
  #filter(grepl("76", Marker)) %>%
  filter(Marker %in% c("pfcrt K76T", "pfcrt K76K/T", "pfcrt 76T", 
                       "pfcrt K76", "pfcrt 76K/T")) %>%
  mutate(year = round((Start.Year + End.Year) / 2, 0),
         # if one or the other is not complete, populate with the value we have:
         year = case_when(is.na(year) & !is.na(Start.Year) ~ Start.Year,
                          is.na(year) & !is.na(End.Year) ~ End.Year,
                          TRUE ~ year)) %>%
  filter(Continent == "Africa") %>%
  suppressWarnings()

plot(crt$Start.Year, crt$Present/crt$Tested, col=as.factor(crt$Marker))
nrow(crt)
sum(crt$Start.Year - crt$End.Year == -1)

# what the heck is happening here
plot(crt$Prevalence...., crt$Prevalence.....incl.mixed)
abline(h=100)

hist(crt$End.Year - crt$Start.Year)

# not at all happy with this but ah well
crt = crt %>%
  # this now only affects one study ... not too worried
  dplyr::summarise(n = dplyr::n(), Present = first(Present),
                   .by = c(uniq_id_publication, Continent, Country, District, Site.Name, Latitude, Longitude,
                           Site.Name.District.Country, Start.Year, End.Year, Tested, PubMedID, Year.Published, Title, Authors, Publication.URL,
                           Contributer, Journal, Estimated.Year, Estimated.Prev, Estimated.Location, Display.on.Surveyor, Mixed.Included, Internal.Notes,
                           Notes, year, Marker)) %>%
  pivot_wider(names_from = Marker, values_from = Present, values_fill = NA) %>%
  # include mixeds if we have them ...
  # what follows is rather silly ...
  filter(!((uniq_id_publication %in% c(614, 1131)  & Mixed.Included == "null") |
             (uniq_id_publication %in% c(1460, 723) & Mixed.Included == "false"))) %>% # these studies were entered twice - once with mixeds included and once without
  mutate(tot = rowSums(across(starts_with("pfcrt")), na.rm = TRUE),
         Present = case_when(!is.na(`pfcrt 76T`) &
                               Mixed.Included == "true" ~ `pfcrt 76T`,
                             !is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "true" ~ `pfcrt 76K/T`,
                             is.na(`pfcrt 76T`) & is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "true" ~ Tested - `pfcrt K76`,
                             !is.na(`pfcrt 76T`) & is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "false" ~ `pfcrt 76T`,
                             !is.na(`pfcrt 76T`) & !is.na(`pfcrt 76K/T`) & 
                               Mixed.Included == "false" ~ `pfcrt 76T` + `pfcrt 76K/T`,
                             is.na(`pfcrt 76T`) & !is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "false" ~ `pfcrt 76K/T`,
                             is.na(`pfcrt 76T`) & is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "false" ~ Tested - `pfcrt K76`,
                             !is.na(`pfcrt 76T`) & !is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "null" ~ `pfcrt 76T` + `pfcrt 76K/T`,
                             !is.na(`pfcrt 76T`) &
                               Mixed.Included == "null" ~ `pfcrt 76T`,
                             is.na(`pfcrt 76T`) & is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "null" ~ Tested - `pfcrt K76`,
                             is.na(`pfcrt 76T`) & !is.na(`pfcrt 76K/T`) &
                               Mixed.Included == "null" ~ `pfcrt 76K/T`,
                             TRUE ~ NA)) %>% # there shouldn't be any NAs left
  filter(Tested > MIN_SAMPLE_SIZE) %>% # check min Tested in k13 script - note in manuscript
  group_by(Longitude, Latitude, year, Tested, uniq_id_publication) %>%
  summarise(n = n(), Present = first(Present)) %>%
  ungroup()

# only one non-conformist ! This is because a mixed has been counted towards all three fields?
which(crt$Present/crt$Tested > 1)
crt$Present[which(crt$Present/crt$Tested > 1)] <- crt$Tested[which(crt$Present/crt$Tested > 1)]

plot(crt$year, crt$Present / crt$Tested)

# check stud size minimum requirement
write.csv(crt, "data/clean/moldm_crt76.csv", row.names = FALSE)



ggplot() + 
  geom_sf(data = filter(world, continent == "Africa"), 
          fill = "white") + 
  geom_point(data = crt %>%
               filter(!is.na(Present)) %>%
               mutate(year_bin = cut_number(year, n = 8)), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_fill_gradientn(colors = iddoPal::iddo_palettes$BlGyRd, 
                       "Prevalence",
                       breaks = c(0, 0.5, 1), 
                       labels = c("0  (all K76)", "0.5", "1  (all 76T)"),
                       limits = c(0,1)) +
  scale_size_continuous(range = c(0.2, 4), trans="sqrt", "Sample size") +
  facet_wrap(~ year_bin, ncol=4) +
  labs(#title = "Pfcrt K76T prevalence") +
    xlab = "Longitude",
    ylab = "Latitude") +
  theme_grey()
ggsave("figures/moldm_crt76.png", width = 4, height = 2.5, scale = 2)


tmp <- crt %>%
  group_by(Longitude, Latitude, year, Tested, uniq_id_publication) %>%
  summarise(n = n(), pick = length(unique(Present)), Present = paste(Present, collapse=",")) %>%
  filter(pick>1)
tmp

# track down what's going on with 179?
# right I'm packing it in - issues with the Mixed.Included value I've put in the filters above
# need to group by lon, lat, year, tested, publication (first(present)) *before* I write to csv

crt %>%
  filter(uniq_id_publication == 1460) %>%
  as.data.frame()

# crt <- mutate(crt,
#               year_bin = cut_number(year, 8))
# 
# ggplot() + 
#   geom_sf(data = st_as_sf(afr), fill = "white") + 
#   #new_scale_color() +
#   geom_point(data = crt, 
#              mapping = aes(x = Longitude, y = Latitude, 
#                            size = Tested,
#                            fill = Present / Tested),
#              col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "Absence", values = c("grey30"), labels=c("")) +
#   scale_fill_viridis_c(name = "Prevalence") +
#   scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
#   facet_wrap(~ year_bin, ncol=4) +
#   labs(title = "Prevalence of PfcrtK76T") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme_grey() 
# ggsave("figs/crt_distn.png", height=3, width=5, scale=2)

