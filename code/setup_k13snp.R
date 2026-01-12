dirs <- list.dirs("output")
snps <- dirs[grepl("k13snp", dirs) & str_count(dirs, "/") == 1]
snps <- str_extract(snps, "(?<=_).*")

mutants <- read.csv("data/clean/moldm_marcse_with_markers.csv")
wts <- read.csv("data/clean/moldm_marcse_wildtypes_to_add.csv")

extract_snp_dataset <- function(snp){
  # filter out matching records of the SNP, create some "absences" where the marker
  # isn't identified (but others are), join with full WT set, write to data/clean/
  
  muts <- filter(mutants, Marker == snp)
  
  no_muts <- mutants %>%
    filter(Marker != snp) %>%
    dplyr::select(Longitude, Latitude, year, Site.Name, Country, Tested) %>% 
    # leaving out tested here:
    anti_join(muts %>% dplyr::select(Longitude, Latitude, year, Site.Name, Country)) %>%
    group_by(Longitude, Latitude, year, Site.Name, Country) %>%
    summarise(Tested = max(Tested)) %>%
    anti_join(wts %>% dplyr::select(Longitude, Latitude, year, Site.Name, Country))
  
  wts <- bind_rows(wts,
                   no_muts) %>%
    mutate(Present = 0)
  
  out <- bind_rows(muts, wts)
  
  print(paste0(snp, "; nrow ", nrow(out), "; npos ", sum(out$Present > 0)))
  
  write.csv(out,
            paste0("data/clean/moldm_marcse_k13snp_", snp, ".csv"))
  
  out
}

suppressMessages(sapply(snps, function(snp){extract_snp_dataset(snp)}))


