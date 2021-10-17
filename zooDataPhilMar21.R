## Claire Davies (CSIRO) 
## zooplankton data for Phil Dyer
## Updated by Phil Dyer

## Created: Mar 2021
## Last Modified: Mar 2021

## suppressPackageStartupMessages({
##   library(lutz)
##   library(lubridate)
##   library(data.table)
##   library(tidyverse)
## })

## rawD <- "RawData"
## outD <- "Output"

#### BGC Phytoplankton #### #################################################################################################################################

# Bring in all BGC zooplankton samples
getBGCSamples <- function(){
BGCSamp <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/BGC_Trip.csv", na = "") %>%
  dplyr::filter(grepl("Z", SAMPLETYPE)) %>%
  dplyr::filter(PROJECTNAME %in% c("NRS", "NRS Ichthyoplankton", "NWS", "IIOE")) %>%
    rename(ProjectName = PROJECTNAME, Station = STATIONNAME, StationCode = STATIONCODE, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateLocal = SAMPLEDATELOCAL,
           BGCcode = TRIP_CODE, Biomass_mgm3 = BIOMASS_MGM3) %>%
  dplyr::filter(!is.na(Latitude) & !is.na(Longitude)) %>%
    mutate(Year = lubridate::year(SampleDateLocal),
           Month = lubridate::month(SampleDateLocal),
           Day = lubridate::day(SampleDateLocal),
           Time_24hr = str_sub(SampleDateLocal, -8, -1),
           tz = tz_lookup_coords(Latitude, Longitude, method = "fast"),
           SampleDateUTC = with_tz(force_tzs(SampleDateLocal, tz, roll = TRUE), "UTC")) 
  return(BGCSamp)
}

getBGCZooData <- function(){
  BGCZdat <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/BGC_Zoop_Raw.csv", na = "") %>%
    rename(BGCcode = TRIP_CODE, TaxonName = TAXON_NAME, Copepod = COPEPOD, TaxonGroup = TAXON_GROUP,
           Genus = GENUS, Species = SPECIES, ZAbund_m3 = ZOOP_ABUNDANCE_M3)
  return(BGCZdat)
}

zoo_process_nrs <- function(){

BGCSamp <- getBGCSamples()
BGCZdat <- getBGCZooData()

# Bring in Change Log
nrsZcl <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/BGC_Zoop_ChangeLog.csv", na = "") %>%
  rename(TaxonName = TAXON_NAME, StartDate = START_DATE, ParentName = PARENT_NAME)

# Check at what level we need change log
nrsclc <- nrsZcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

BGCCop1 <- BGCZdat %>%
  filter(!TaxonName %in% levels(as.factor(nrsclc$TaxonName)) &  (TaxonGroup %in% c('Copepod', 'Ciliophora', 'Cladoceran') | grepl("Dolio", TaxonName)) &
         Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(ProjectNumber = 599, 
         TaxonGroup = ifelse(TaxonGroup == 'Ciliophora', 'Tintinnid', TaxonGroup),
         Species = paste0(word(Genus,1)," ", word(Species,1))) %>% # bin complexes
  group_by(BGCcode, TaxonGroup, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop")

BGCCop1 <- left_join(BGCCop1, BGCSamp, by = "BGCcode") %>%
  mutate(StartDate = ymd("2007-12-19")) %>%  # avoids nulls in pivot
  group_by(BGCcode, Station, Latitude, Longitude, SampleDateUTC, TaxonGroup, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
BGCCop2 <- BGCZdat %>%
  filter(TaxonName %in% levels(as.factor(nrsclc$TaxonName)) &  (TaxonGroup %in% c('Copepod', 'Ciliophora', 'Cladoceran') | grepl("Dolio", TaxonName)) &
           Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species)) %>% 
  mutate(ProjectNumber = 599, 
         TaxonGroup = ifelse(TaxonGroup == 'Ciliophora', 'Tintinnid', TaxonGroup),
         Species = paste0(word(Genus,1), " ", word(Species,1))) %>% # bin complexes
  left_join(nrsZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% 
  drop_na(Species) %>%
  group_by(BGCcode, StartDate, TaxonGroup, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(BGCCop2$Species)) {
  Taxon <- BGCCop2 %>% dplyr::select(Species) %>% unique()
  Taxon <- as.character(Taxon$Species[i] %>% droplevels())
  
  Dates <- as.data.frame(BGCCop2) %>%
    filter(Species == Taxon) %>% 
    slice(1) %>% 
    droplevels()
  
  copes1 <- as.data.frame(BGCCop2) %>%
    filter(Species == Taxon) %>% 
    droplevels() 
  
  copes <- left_join(BGCSamp, copes1, by = "BGCcode") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           TaxonGroup = replace(TaxonGroup, is.na(TaxonGroup), unique(copes1$TaxonGroup)),
           ZAbund_m3 = replace(ZAbund_m3, is.na(ZAbund_m3) & StartDate>SampleDateLocal, -999)) %>% 
    drop_na(ZAbund_m3) %>%
    group_by(BGCcode, Station, Latitude, Longitude, SampleDateUTC, TaxonGroup, Species) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
    as.data.frame()     
  BGCCop1 <- rbind(BGCCop1, copes)
}

BGCRawZ1 <- BGCCop1 %>%
  group_by(Latitude, Longitude, SampleDateUTC, TaxonGroup, Species) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3), .groups = "drop") %>% 
  mutate(ProjectNumber = 599,
         TaxonGroup = ifelse(TaxonGroup == 'Ciliophora', 'Tintinnid', TaxonGroup)) %>%
  arrange(-desc(Species)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified

}
zoo_process_cpr <- function(){
# Bring in CPR zoo
# Bring in all CPR zooplankton samples
cprZsamp <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/CPR_Samp.csv", na = "") %>%
  dplyr::filter( grepl("Z", SAMPLETYPE)) %>%
  rename(Sample = SAMPLE, BGCcode = TRIP_CODE, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLEDATEUTC) %>%
  mutate(Year = lubridate::year(SampleDateUTC),
         Month = lubridate::month(SampleDateUTC),
         Day = lubridate::day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1)) # hms doesn"t seem to work on 00:00:00 times

# Bring in plankton summary data
cprZdat <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/CPR_Zoop_Raw.csv", na = "") %>%
  rename(Sample = SAMPLE, TaxonName = TAXON_NAME, Copepod = TAXON_GROUP, TaxonGroup = TAXON_GRP01,
         Genus = GENUS, Species = SPECIES, ZAbund_m3 = ZOOP_ABUNDANCE_M3)

# Bring in Change Log
cprZcl <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/CPR_Zoop_ChangeLog.csv", na = "") %>%
  rename(TaxonName = TAXON_NAME, StartDate = STARTDATE, ParentName = PARENT_NAME)

# Check at what level we need change log
clc <- cprZcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

cprCop1 <- cprZdat %>% 
  filter(!TaxonName %in% levels(as.factor(clc$TaxonName)) & Species != "spp." & !grepl("cf.", Species) & !grepl("grp", Species) & !is.na(Species) &
           (TaxonGroup %in% c('Copepod', 'Ciliophora', 'Cladoceran') | grepl("Dolio", TaxonName))) %>%
  mutate(Species = paste0(word(Genus,1)," ", word(Species,1))) %>% # bin complexes
  group_by(Sample, Species, TaxonGroup) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop")

cprCop1 <- left_join(cprCop1, cprZsamp, by = "Sample") %>%
  drop_na(TaxonGroup) %>%
  mutate(StartDate = ymd("2007-12-19")) %>% # avoids nulls in pivot
  group_by(Latitude, Longitude, SampleDateUTC, TaxonGroup, Species) %>%
  summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprCop2 <- cprZdat %>% 
  filter(TaxonName %in% levels(as.factor(clc$TaxonName)) 
         & Species != "spp." & !is.na(Species) & !grepl("cf.", Species) & !grepl("grp", Species) &
           (TaxonGroup %in% c('Copepod', 'Ciliophora', 'Cladoceran') | grepl("Dolio", TaxonName))) %>% 
  mutate(Species = paste0(word(Genus,1)," ", word(Species,1))) %>% # bin complexes
  left_join(cprZcl, by = "TaxonName") %>%
  mutate(Species = as_factor(Species)) %>% drop_na(Species) %>%
  group_by(Sample, StartDate, Species, TaxonGroup) %>% 
  summarise(ZAbund_m3 = sum(ZAbund_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprCop2$Species)) {
  Spe <- cprCop2 %>% dplyr::select(Species) %>% unique()
  Spe <- as.character(Spe$Species[i] %>% droplevels())
  
  Dates <- as.data.frame(cprCop2) %>% 
    filter(Species == Spe) %>% 
    slice(1) %>% 
    droplevels()
  
  copes1 <- as.data.frame(cprCop2) %>% 
    filter(Species == Spe) %>%
    droplevels() 
  
  copes <- cprZsamp %>% left_join(copes1, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           Species = replace(Species, is.na(Species), Dates$Species),
           TaxonGroup = replace(TaxonGroup, is.na(TaxonGroup), unique(copes1$TaxonGroup)),
           ZAbund_m3 = replace(ZAbund_m3, is.na(ZAbund_m3) & StartDate>SampleDateUTC, -999)) %>% 
    drop_na(ZAbund_m3) %>%
    group_by(Latitude, Longitude, SampleDateUTC, TaxonGroup, Species) %>%
    summarise(ZAbund_m3 = sum(ZAbund_m3), .groups = "drop") %>% 
    as.data.frame()     
  cprCop1 <- rbind(cprCop1, copes)
}

cprRawZ1 <- cprCop1 %>% 
  group_by(Latitude, Longitude, SampleDateUTC, TaxonGroup, Species) %>%
  summarise(ZAbund_m3 = max(ZAbund_m3), .groups = "drop") %>% 
  arrange(-desc(Species)) %>% 
  mutate(ProjectNumber = 597,
         TaxonGroup = ifelse(TaxonGroup == 'Ciliophora', 'Tintinnid', TaxonGroup)) %>%
  as.data.frame() 

CPRRawZ2 <- left_join(cprZsamp, cprZdat, by = "Sample") %>% #samples with no taxa associated with them
  filter(is.na(TaxonGroup)) %>%
  dplyr::select(c(Latitude, Longitude, SampleDateUTC, TaxonGroup)) %>% unique() %>%
  mutate(Species = 'No taxa', 
         TaxonGroup = 'No taxa',
         ZAbund_m3 = NA, 
         ProjectNumber = 597)

cpr_all <- rbind(cprRawZ1, CPRRawZ2)
}

zoo_process_other <- function(root_dir = "./"){
# Bring in other data
otherZ <- read_csv(glue::glue("{root_dir}/zoop_phil_feb21.csv"), na = "") %>%
  rename(Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLE_TIME_UTC, TaxonName = TAXON_NAME, 
         ZAbund_m3 = ABUNDANCE_M3, ProjectNumber = PROJECTNUMBER, TaxonGroup = FUNCTIONAL_GROUP) %>%
  dplyr::select(Latitude , Longitude, SampleDateUTC , TaxonGroup, TaxonName, ZAbund_m3 , ProjectNumber) %>%
  filter(!grepl("invalid", TaxonName) &!grepl("complex", TaxonName)) %>%
  mutate(TaxonGroup = ifelse(TaxonGroup %in% c('Calanoid', 'Cyclopoid', 'Poecilostomatoida', 'Harpacticoid'), 'Copepod', TaxonGroup), 
         Species = gsub("\\s*\\([^\\)]+\\)","", TaxonName)) %>%
  dplyr::select(-TaxonName)
}

load_zoo_data <- function(root_dir = "./") {
zoop <- rbind(zoo_process_nrs(), zoo_process_cpr(), zoo_process_other(root_dir))
}
