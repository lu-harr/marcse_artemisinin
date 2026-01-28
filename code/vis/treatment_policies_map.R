library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(tidyverse)

# country shps for plotting/masking
world <- ne_countries(scale="medium", returnclass = "sf")

afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

# national first-line treatment policies from 2025 world malaria report
dat <- read_xlsx("world_malaria_report_treatment_policies.xlsx")
dat$country
dat$`Uncomplicated confirmed`

confirmed <- dat %>%
  separate_longer_delim(`Uncomplicated confirmed`, delim = ";") %>%
  mutate(uncomp_conf = trimws(`Uncomplicated confirmed`),
         uncomp_conf = ifelse(uncomp_conf == "NA", NA, uncomp_conf),
         uncomp_conf = ifelse(uncomp_conf == "AL-PQ", "AL+PQ", uncomp_conf)) %>%
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from = uncomp_conf,
    values_from = value,
    values_fill = FALSE
  ) %>%
  mutate(country = case_when(country == "Central African Republic" ~ "Central African Rep.",
                             country == "Democratic Republic of the Congo" ~ "Dem. Rep. Congo",
                             country == "Equatorial Guinea" ~ "Eq. Guinea",
                             country == "Eswatini" ~ "eSwatini",
                             country == "Sao Tome and Principe" ~ "São Tomé and Principe",
                             country == "South Sudan" ~ "S. Sudan",
                             country == "Mainland" ~ "Tanzania",
                             TRUE ~ country)) %>%
  mutate(al_asaq = case_when(AL & !`AS+AQ` ~ "AL",
                             !AL & `AS+AQ` ~ "AS+AQ",
                             NA ~ "")) # not doing the others as this is just for legend
         
confirmed <- right_join(afr, confirmed, by = join_by(name == country))

# tanzania is a problem to map .. diff policies for mainland/zanzibar

pal = c("yellow", "#43a2ca", "#f03b20")

library(ggpattern)

ggplot() +
  # geom_sf(data = afr) + # include to show borders of all countries
  geom_sf(data = confirmed %>% filter(!is.na(al_asaq)),
          aes(fill = al_asaq)) +
  geom_sf_pattern(data = confirmed %>% 
                    filter(`DHA-PPQ` | `AL+PQ`) %>%
                    mutate(flag = case_when(`DHA-PPQ` ~ "DHA-PPQ",
                                            `AL+PQ` ~ "AL+PQ",
                                            TRUE ~ NA)),
                  aes(pattern = flag),
                  fill = NA,
                  pattern_fill = "grey50",
                  pattern_spacing = 0.02,
                  pattern_density = 0.001) +
  scale_pattern_manual(name = "", 
                       values = c("DHA-PPQ" = "circle", 
                                  "AL+PQ" = "stripe")) +
  scale_fill_manual(name = "", 
                    values = c("AL" = alpha(pal[1], 0.3),
                               "AS+AQ" = alpha(pal[2], 0.3))) +
  labs(title = "First-line treatment policies (2025)")
ggsave("~/Desktop/test.png")

# needs also: AS-PY, AS+MQ, etc .........






