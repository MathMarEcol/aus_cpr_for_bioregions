## Claire Davies (CSIRO) 
## phytoplankton data for Phil Dyer
## Updated by Phil Dyer

## Created: Feb 2021
## Last Modified: Mar 2021

## suppressPackageStartupMessages({
##   library(lutz)
##   library(lubridate)
##   library(tidyverse)
## })

## rawD <- "RawData"
## outD <- "Output"

#### BGC Phytoplankton #### #################################################################################################################################
phyto_process_nrs <- function(){
# Bring in all BGC phytoplankton samples
BGCPsamp <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/c6cbc0166a4d1c13f336923bac3dbe052eaad172/Plankton/RawData/BGC_Trip.csv", na = "") %>%
  dplyr::filter( grepl("P", SAMPLETYPE)) %>%
  dplyr::filter(PROJECTNAME %in% c("NRS", "NRS Ichthyoplankton", "NWS", "IIOE")) %>%
    rename(ProjectName = PROJECTNAME, Station = STATIONNAME, StationCode = STATIONCODE, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateLocal = SAMPLEDATELOCAL,
           BGCcode = TRIP_CODE, Biomass_mgm3 = BIOMASS_MGM3) %>%
  dplyr::filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  mutate(tz = tz_lookup_coords(Latitude, Longitude, method = "fast"),
          SampleDateUTC = with_tz(force_tzs(SampleDateLocal, tz, roll = TRUE), "UTC")) 

# Bring in plankton data
BGCPdat <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/c6cbc0166a4d1c13f336923bac3dbe052eaad172/Plankton/RawData/BGC_Phyto_Raw.csv", na = "") %>%
  rename(BGCcode = TRIP_CODE, TaxonName = TAXON_NAME, TaxonGroup = TAXON_GROUP, Genus = GENUS, Species = SPECIES,
         Cells_L = CELL_L, Biovolume_uM3_L = BIOVOLUME_UM3L)

# Bring in Change Log
BGCPcl <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/c6cbc0166a4d1c13f336923bac3dbe052eaad172/Plankton/RawData/BGC_Phyto_ChangeLog.csv", na = "", col_select = 1:3) %>%
  rename(TaxonName = TAXON_NAME, StartDate = STARTDATE, ParentName = PARENT_NAME)

#### Species Abund ####

# Check at what level we need change log
nrsls <- BGCPcl %>%
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species

BGCSpecP1 <- BGCPdat %>%
  filter(!TaxonName %in% levels(as.factor(nrsls$TaxonName)) &
         !TaxonGroup %in% c("Ciliate", "Other") &
         Species != "spp." & !is.na(Species) & !grepl("/", Species) & !grepl("cf.", Species) & !grepl("complex", Species) & !grepl("type", Species)) %>%
  mutate(TaxonName = paste0(word(Genus,1), ' ', word(Species,1))) %>% 
  group_by(BGCcode, TaxonName) %>%
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop")

BGCSpecP1 <-  left_join(BGCSpecP1, BGCPsamp, by = "BGCcode") %>%
  mutate(StartDate = ymd("2007-12-19"),
         Cells_L = ifelse(is.na(Cells_L), 0, Cells_L)) %>% 
  group_by(BGCcode, Station, Latitude, Longitude, SampleDateUTC, TaxonName) %>%
  summarise(Cells_L = sum(Cells_L), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
BGCSpecP2 <- BGCPdat %>%
  filter(TaxonName %in% levels(as.factor(nrsls$TaxonName)) &
         !TaxonGroup %in% c("Ciliate", "Other") &
         Species != "spp." & !is.na(Species) & !grepl("/", Species) & !grepl("complex", Species) & !grepl("type", Species)) %>% 
  left_join(BGCPcl, by = "TaxonName") %>%
  mutate(TaxonName = paste0(word(Genus,1), ' ', word(Species,1)),
         TaxonName = as_factor(TaxonName)) %>% 
  drop_na(TaxonName) %>%
  group_by(BGCcode, StartDate, TaxonName) %>%
  summarise(Cells_L = sum(Cells_L, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(BGCSpecP2$TaxonName)) {
  Taxon <- BGCSpecP2 %>% dplyr::select(TaxonName) %>% unique()
  Taxon <- as.character(Taxon$TaxonName[i] %>% droplevels())
  
  Dates <- as.data.frame(BGCSpecP2) %>%
    filter(TaxonName == Taxon) %>% 
    slice(1) %>% 
    droplevels()
  
  spec1 <- as.data.frame(BGCSpecP2) %>%
    filter(TaxonName == Taxon) %>% 
    droplevels() 

  spec <- BGCPsamp %>%
    left_join(spec1, by = "BGCcode") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           Cells_L = replace(Cells_L, is.na(Cells_L) & StartDate>SampleDateLocal, -999)) %>% 
    drop_na(Cells_L) %>%
    group_by(BGCcode, Station, Latitude, Longitude, SampleDateUTC, TaxonName) %>%
    summarise(Cells_L = sum(Cells_L), .groups = "drop") %>% 
    as.data.frame()     
  BGCSpecP1 <- rbind(BGCSpecP1, spec)
}

BGCPhyto <- BGCSpecP1 %>%
  group_by(BGCcode, Station, Latitude, Longitude, SampleDateUTC, TaxonName) %>%
  summarise(Cells_L = max(Cells_L), .groups = "drop") %>% 
  dplyr::select(Latitude, Longitude, SampleDateUTC, TaxonName, Cells_L) %>%
  mutate(ProjectNumber = 599) %>%
  arrange(-desc(TaxonName)) %>% 
  as.data.frame() 
# select maximum value of duplicates, but leave -999 for all other occurences as not regularly identified
}

#### CPR Phytoplankton #######################################################################################################################################################
phyto_process_cpr <- function(){
# Bring in all CPR phytoplankton samples
cprPsamp <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/c6cbc0166a4d1c13f336923bac3dbe052eaad172/Plankton/RawData/CPR_Samp.csv", na = "") %>%
  dplyr::filter( grepl("P", SAMPLETYPE)) %>%
## cprPsamp <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/9b4f5865af778164fbefc480913bcca3b17be127/Plankton/RawData/PsampCPR.csv", na = "(null)") %>%
  rename(BGCcode = TRIP_CODE, Sample = SAMPLE, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLEDATEUTC) %>%
  mutate(Year = year(SampleDateUTC),
         Month = month(SampleDateUTC),
         Day = day(SampleDateUTC),
         Time_24hr = str_sub(SampleDateUTC, -8, -1)) # hms doesn"t seem to work on 00:00:00 times

# Bring in plankton data
cprPdat <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/c6cbc0166a4d1c13f336923bac3dbe052eaad172/Plankton/RawData/CPR_Phyto_Raw.csv", na = "") %>%
  rename(Sample = SAMPLE, TaxonName = TAXON_NAME, TaxonGroup = TAXON_GROUP, Genus = GENUS, Species = SPECIES, PAbun_m3 = PHYTO_ABUNDANCE_M3, BioVolume_um3_m3 = BIOVOL_UM3M3)

# Bring in Change Log
cprPcl <- read_csv("https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/c6cbc0166a4d1c13f336923bac3dbe052eaad172/Plankton/RawData/CPR_Phyto_ChangeLog.csv", na = "") %>%
  rename(TaxonName = TAXON_NAME, StartDate = STARTDATE, ParentName = PARENT_NAME)

#### CPR PHYTO ABUND SPECIES ####

# Check at what level we need change log
cls <- cprPcl %>% 
  mutate(same = ifelse(TaxonName == ParentName, "yes", "no")) %>%
  filter(same == "no") # no changes at genera level

# for non change log species
cprSpecP1 <- cprPdat %>% 
  filter(!TaxonName %in% levels(as.factor(cls$TaxonName)) &
           !TaxonGroup %in% c("Ciliate", "Other") &  
         Species != "spp." & !is.na(Species) & !grepl("/", Species)  & !grepl("cf.", Species) & !grepl("complex", Species) & !grepl("type", Species)) %>% 
  mutate(TaxonName = paste0(word(Genus,1), ' ', word(Species,1))) %>% 
  group_by(Sample, TaxonName) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop")

cprSpecP1 <- left_join(cprSpecP1, cprPsamp, by = "Sample") %>%
  mutate(StartDate = ymd("2007-12-19"))  %>% 
  group_by(BGCcode, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(PAbun_m3 = sum(PAbun_m3), .groups = "drop") %>% 
  as.data.frame()

# add change log species with -999 for NA"s and real absences as 0"s
cprSpecP2 <- cprPdat %>% 
  filter(TaxonName %in% levels(as.factor(cls$TaxonName)) &
           !TaxonGroup %in% c("Ciliate", "Other") &
         Species != "spp." & !is.na(Species)  & !grepl("cf.", Species) & !grepl("complex", Species) & !grepl("type", Species)
         & !grepl("/", Species)) %>%
  mutate(TaxonName = paste0(word(Genus,1), ' ', word(Species,1))) %>% 
  left_join(cprPcl, by = "TaxonName") %>%
  mutate(TaxonName = as_factor(TaxonName)) %>% 
  drop_na(TaxonName) %>%
  group_by(Sample, StartDate, TaxonName) %>% 
  summarise(PAbun_m3 = sum(PAbun_m3, na.rm = TRUE), .groups = "drop") 

for (i in 1:nlevels(cprSpecP2$TaxonName)) {
  Taxon <- cprSpecP2 %>% dplyr::select(TaxonName) %>% unique()
  Taxon <- as.character(Taxon$TaxonName[i] %>% droplevels())

    Dates <- as.data.frame(cprSpecP2) %>% 
    filter(TaxonName == Taxon) %>% 
    slice(1) %>% 
    droplevels()
  
  spec <- as.data.frame(cprSpecP2) %>% 
    filter(TaxonName == Taxon) %>% 
    droplevels() 
  
  spec <- cprPsamp %>% 
    left_join(spec, by = "Sample") %>%
    mutate(StartDate = replace(StartDate, is.na(StartDate), Dates$StartDate),
           TaxonName = replace(TaxonName, is.na(TaxonName), Dates$TaxonName),
           PAbun_m3 = replace(PAbun_m3, is.na(PAbun_m3) & StartDate>SampleDateUTC, -999)) %>% 
    group_by(BGCcode, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
    drop_na(PAbun_m3) %>%
    summarise(PAbun_m3 = sum(PAbun_m3), .groups = "drop") %>% 
    as.data.frame()     
  cprSpecP1 <- rbind(cprSpecP1, spec)
}

CPRPhyto <- cprSpecP1 %>% 
  group_by(BGCcode, Latitude, Longitude, SampleDateUTC, Year, Month, Day, Time_24hr, TaxonName) %>%
  summarise(Cells_L = max(PAbun_m3), .groups = "drop") %>% 
  dplyr::select(Latitude, Longitude, SampleDateUTC, TaxonName, Cells_L) %>%
  mutate(Cells_L = ifelse(Cells_L == -999, -999, Cells_L/1000),
         ProjectNumber = 597) %>%
  arrange(-desc(TaxonName)) %>% 
  as.data.frame() 

CPRRawP2 <- left_join(cprPsamp, cprPdat, by = "Sample") %>% #samples with no taxa associated with them
  dplyr::select(c(Latitude, Longitude, SampleDateUTC, TaxonName, PAbun_m3)) %>%
  rename(Cells_L = PAbun_m3) %>%
  filter(is.na(TaxonName)) %>%
  mutate(TaxonName = 'No taxa', 
         ProjectNumber = 597)

CPRPhyto <- rbind(CPRPhyto, CPRRawP2)
}
## Steve Brett Data
phyto_process_brett <- function(root_dir = "./"){
  ## Available at https://raw.githubusercontent.com/MathMarEcol/HABS_STAR/9b4f5865af778164fbefc480913bcca3b17be127/Brett_data_feb2020.csv,
  ## but the repo is private, so I am providing the file locally.
BrettDat <- read_csv(glue::glue("{root_dir}/Brett_data_feb2020.csv"), na = c("NA",""),
                     col_types = cols(LATITUDE = col_double(),
                                      LONGITUDE = col_double())) %>%
  rename(Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLE_DATE, TaxonName = TAXON_NAME, Cells_L = ABUNDANCE) %>%
  filter(!is.na(Latitude) & !is.na(Longitude) & !grepl(" cyst", TaxonName) &
         !grepl("spp.", TaxonName) & !grepl("/", TaxonName)  & !grepl("cf.", TaxonName) & !grepl("complex", TaxonName) & !grepl("type", TaxonName) & !grepl("iatom", TaxonName) &
         !grepl("lagellate", TaxonName)& !grepl("uglenida", TaxonName) &!grepl("yanophyt", TaxonName)) %>%
  dplyr::select(Latitude, Longitude, SampleDateUTC, TaxonName, Cells_L) %>%
  mutate(TaxonName = ifelse(grepl("width", TaxonName), paste0(word(TaxonName, 1)," ", word(TaxonName, 2)), TaxonName),
         ProjectNumber = 794, SampleDateUTC = as.POSIXct(lubridate::as_date( SampleDateUTC)))

}

## Other Data 
phyto_process_other <- function(root_dir = "./"){
OtherDat <- read_csv(glue::glue("{root_dir}/phyto_other_phil_feb21.csv"), na = "") %>%
  dplyr::rename(Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateUTC = SAMPLE_TIME_UTC, TaxonName = TAXON_NAME, Cells_L = CELLS_L,
         ProjectNumber = PROJECT_ID) %>%
  dplyr::filter(!grepl("/", TaxonName)& !grepl("spp.", TaxonName)& !grepl("cf.", TaxonName) & !grepl(" cyst", TaxonName)) %>%
  dplyr::select(Latitude, Longitude, SampleDateUTC, TaxonName, Cells_L, ProjectNumber) %>%
  dplyr::mutate(TaxonName = ifelse(grepl("Richelia", TaxonName), paste0(word(TaxonName, 1)," ", word(TaxonName, 2)), TaxonName)) %>%
  dplyr::mutate(TaxonName = gsub("\\s*\\([^\\)]+\\)","", TaxonName))
}


load_phyto_data <- function(root_dir = "./"){

  bound <- rbind(phyto_process_other(root_dir),
phyto_process_brett(root_dir),
phyto_process_cpr(),
phyto_process_nrs()
  )
}
