# read in moldm and clean up into with_markers and aggregated formats

# this script does all of our packages and brings in the WHO markers list:
source("code/data_clean/format_moldm.R")

moldm <- format_moldm_k13("data/raw/db_20260211/novartis.csv") %>%
  # do a quick tidy up of PubMedIDs
  clean_up_pmids()
  
# what's left?
moldm %>%
filter(str_length(PubMedID) != 8 & PubMedID != "Unpublished") %>%
dplyr::select(PubMedID, Title) %>%
arrange(Title) %>%
unique() %>%
as.vector()
# should be a doi or "Unpublished"

# notes on double mutants:
# 26667053: double mutants not counted in single mutant counts
# 28039354: double mutants not counted in single mutant counts
# 39136465: double mutants not counted in single mutant counts
# 40439506: double mutants not counted in single mutant counts
# 40790052: double mutants not counted in single mutant counts



write.csv(moldm %>%
            filter(Present / Tested <= 1 & Tested > MIN_SAMPLE_SIZE), 
          # the Tested filter should make a difference but the Prevalence shouldn't
          "data/clean/moldm_with_markers.csv")

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

plot(mutants$year, mutants$Present/mutants$Tested, xlab="Year", ylab="Prevalence")

wildtypes <- moldm %>%
  filter(Marker == "k13 wildtype") %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(unique(Present)), n=n(), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>% 
  # check how many simultaneous wildtype entries we have?
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country, Present, pubs) %>%
  ungroup() %>%
  suppressMessages()

# possible duplicates?
wildtypes[which(wildtypes$Present/wildtypes$Tested > 1),
          c("Longitude", "Latitude",  "year", "Tested", "Site.Name",
            "Present", "pubs")]
# it doesn't matter so much because of how we add them below:

wildtypes_to_add <- anti_join(
  # locations in `wildtypes` that do not occur in `mutants`,
  # paying attention to `Tested` but NOT to `pubs`
  wildtypes %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country),
  mutants %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country)) %>%
  mutate(Present = 0) %>%
  suppressMessages()

message(paste("Number of rows of wildtypes to add:", nrow(wildtypes_to_add)))

with_wildtypes <- full_join(mutants, wildtypes_to_add) %>%
  filter(Tested > MIN_SAMPLE_SIZE) %>%
  suppressMessages()

write.csv(with_wildtypes %>%
            dplyr::select(-c(pubs)),
          "data/clean/moldm_k13_nomarker.csv",
          row.names = FALSE)

message(paste("Number of testees:", 
              with_wildtypes %>% 
                group_by(year, Longitude, Latitude, Tested) %>% 
                summarise(n=n()) %>% 
                ungroup() %>%
                dplyr::select(Tested) %>% 
                sum()))

###############################################################################
# all vis down here
with_wildtypes <- read.csv("data/clean/moldm_k13_nomarker.csv")
moldm <- read.csv("data/clean/moldm_with_markers.csv")

plot(with_wildtypes$year, with_wildtypes$Present/with_wildtypes$Tested, 
     xlab="Year", ylab="Prevalence")
ggplot(data = moldm %>% filter(mutant)) +
  geom_point(aes(x = year, y = Present / Tested, size = Tested))
moldm %>% filter(Present/Tested > 0.4 & mutant)


to_vis <- with_wildtypes %>% 
  mutate(year_bin = cut(year, breaks = c(min(year)-1, 2009, 2012, 2015, 2018, 2021, max(year)))) %>%
  # mutate(year_bin = case_when(year <= 2007 ~ "1995-2009",
  #                             year <= 2011 ~ "2010-2011",
  #                             year <= 2013 ~ "2012-2013",
  #                             year <= 2015 ~ "2014-2015",
  #                             year <= 2017 ~ "2016-2017",
  #                             year <= 2019 ~ "2018-2019",
  #                             year <= 2021 ~ "2020-2021",
  #                             TRUE ~ "2022-2024")) %>%
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
  scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
  facet_wrap(~ year_bin, ncol=3) +
  labs(title = "Prevalence of k13 markers in Africa") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_grey() 
# ggsave("figures/archive/abs_grey.png", height=3, width=5, scale=2)

markers <- moldm %>%
  filter(mutant) %>% # no filter on Tested
  group_by(year, strip_marker) %>%
  summarise(present = sum(Present)) %>%
  rename(marker = strip_marker)

# put test intensity in the background
bg <- moldm %>%
  group_by(Longitude, Latitude, year, PubMedID) %>%
  summarise(Tested = sum(unique(Tested))) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(Tested = sum(Tested))

# I cannot state again how intensely easy this is in base graphics

markers_keep <- markers %>%
  filter(present > 2) %>%
  arrange(present) %>%
  ungroup() %>%
  dplyr::select(marker) %>%
  unique()

npal = 6
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
  bind_rows(wts %>% mutate(Longitude = as.numeric(Longitude),
                           Latitude = as.numeric(Latitude))) %>% #  %>% mutate(Marker = "wt")
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

bg_scale <- 50
bg_col <- "grey65"

library(viridisLite)
 p1 <- 
ggplot() +
  geom_bar(data = bg, aes(x = year, y = Tested / bg_scale), 
           stat = "identity", fill = "grey75") +
  geom_line(data = df, size = 0.7,
            aes(x = year, y = present, group = marker, color = marker, linetype = marker)) +
  geom_point(data = df, 
             aes(x = year, y = present, group = marker, color = marker)) +
  labs(
    title = "(a) Detected mutations by year",
    color = "Marker",
    linetype = "Marker"
  ) +
  scale_color_manual(values = rep(c(viridis(5), "#E37210", iddoblue, "#c7047c"), 2)) +
  scale_linetype_manual(values = rep(1:2, each = 8)) +
  scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests")) +
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
  scale_x_continuous(breaks = 2005:2025) 
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
             col = "grey50", pch=21, stroke = 0.2) +
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
  theme(plot.background = element_rect(fill='transparent', color=NA))

library(grid)

rect <- rectGrob(
  x = 0.915,
  y = 0.27,
  width = unit(0.2, "in"),
  height = unit(0.2, "in"),
  gp = gpar(fill = "grey70", colour="grey50", alpha = 0.5)
)

p2 <- ggdraw(p2) +
  draw_grob(rect) + 
  draw_label("Absence", x = 0.96, y = 0.27, size = 10)
#theme(#strip.background = element_blank(),
#strip.text.x = element_blank(),
#strip.placement = "outside") +
#ggsave("figs/markers_disagg.png", height = 4, width = 5, scale = 2)

plot_grid(p1, p2, ncol = 1, rel_heights = c(0.6, 2))

# p1 + p2 + plot_layout(ncol = 1, heights = c(0.5,2))

# would be sick if I could make the word "absence" bigger but I've had enough of 
# snarky people on ggplot stack overflow for about three years
ggsave("figures/markers_disagg.png", height = 5, width = 5.1, scale = 2)




