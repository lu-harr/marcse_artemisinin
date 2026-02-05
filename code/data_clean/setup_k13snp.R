# do this after you have run merge MARCSE script

dirs <- list.dirs("output")
snps <- dirs[grepl("k13snp", dirs) & str_count(dirs, "/") == 1]
snps <- str_extract(snps, "(?<=_).*")

mutants <- read.csv("data/clean/moldm_marcse_with_markers.csv")
wts <- read.csv("data/clean/moldm_marcse_wildtypes_to_add.csv")

extract_snp_dataset <- function(snp){
  # filter out matching records of the SNP, join with full WT set, write to data/clean/
  
  muts <- filter(mutants, Marker == snp)
  
  wts <- anti_join(wts, 
                   muts %>% dplyr::select(Longitude, Latitude, year, 
                                          Site.Name, Country)) %>%
    group_by(Longitude, Latitude, year, Site.Name, Country) %>%
    summarise(Tested = max(Tested)) %>%
    mutate(Present = 0)
  
  out <- bind_rows(muts, wts)
  
  print(paste0(snp, "; nrow ", nrow(out), "; npos ", sum(out$Present > 0)))
  
  write.csv(out,
            paste0("data/clean/moldm_marcse_k13snp_", snp, ".csv"))
}

suppressMessages(sapply(snps, function(snp){extract_snp_dataset(snp)}))


# having a look at Southern cluster
ggplot(data = mutants %>%
         filter(Present > 0 & !is.na(status)) %>%
         mutate(year_bin = cut(year, breaks = c(min(year) - 1, 2012, 2016, 2020, 2025)))) +
  geom_sf(data = afr) +
  geom_point(aes(x = Longitude, y = Latitude, #col = PubMedID, 
                 size = Present/Tested), alpha = 0.2) +
  facet_wrap(~year_bin)

ggplot(data = mutants %>%
         filter(Marker == "A675V" & Present > 0 & year < 2018) %>%
         mutate(year_bin = cut(year, 
                               breaks = c(min(year) -1, 2016, 2018, 2020, 2022, 2025)))) +
  geom_sf(data = afr) +
  geom_point(aes(x = Longitude, y = Latitude, col = PubMedID, 
                 size = Tested), alpha = 0.2) +
  facet_wrap(~year_bin)
ggsave("~/Desktop/test.png", height = 10, width = 15)


ggplot(data = mutants %>%
         filter(Present > 0 & Latitude < -10 & !is.na(status) & year > 2022)) + #%>%
         #mutate(year_bin = cut(year, breaks = c(min(year) - 1, 2022, 2023, 2025)))) +
  geom_sf(data = afr) +
  geom_point(aes(x = Longitude, y = Latitude, col = Marker, 
                 size = Tested), alpha = 0.2) +
  facet_wrap(~Marker, ncol = 2)
ggsave("~/Desktop/test.png", height = 10, width = 10)


mutants %>%
  filter(Marker == "A675V" & Present > 0 & 
           year > 2012 & year <= 2016) %>%
           #year > 2016) %>%
  dplyr::select(Title, PubMedID) %>%
  unique()
# leaving out Wang et al - imported cases into China

mutants %>%
  filter(Marker == "A675V" & Present > 0 & year < 2018) %>%
  dplyr::select(Title, PubMedID) %>%
  unique()

mutants %>%
  filter(Marker == "P441L" & Present > 0) %>%
  filter(PubMedID %in% c("40744006", "40744004"))

mutants %>%
  filter(Marker == "P574L" & Present > 0) %>%
  filter(PubMedID %in% c("40744006", "40744004"))

mutants %>%
  filter(Country == "Nigeria" & Present > 0 & Marker != "wildtype" & !is.na(status))

