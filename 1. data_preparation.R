##### Ecoregion - Data preparation ####
library(tidyverse)
library(sf)
library(raster)

## GIS DATA ####
# ECOREGIONS - load only palearctic ecoregions (to reduce time and memory use)
ecoregions <- st_read('data_raw/Ecoregions2017/Ecoregions2017.shp',
                      query="SELECT * FROM Ecoregions2017 WHERE REALM = 'Palearctic'",
                      stringsAsFactors=FALSE)

# Now switch to a Lambert Azimuthal Equal Area projection as recommended by EEA
# http://epsg.io/3035 for area calculations and analysis
ecoregions <- st_transform(ecoregions, crs=3035)
ecoregions$area <- st_area(ecoregions)
ecoregions$area_sqkm <- as.numeric(ecoregions$area) / 1e6

# LAND MAP 
land <- st_read('data_raw/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')
land <- st_transform(land, crs=3035)

# FOCAL AREA
# create a region for the study area				  
study_bbox <- st_bbox(c(xmin=2500000, xmax=5900000, ymin=1500000, ymax=5400000), crs=3035)
study <- st_as_sfc(study_bbox)

# Can now reduce the ecoregions back to the ones that intersect the study area
ecoregions <- ecoregions[st_intersects(ecoregions, study, sparse=FALSE),]
land <- land[st_intersects(land, study, sparse=FALSE),]

## STUDY DATA ----
# SITES
sites <- read.csv(file = "data_raw/sites.csv", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

## Add climate data 
# load the data
MAT <- stack(raster("data_raw/wc_MAT_30s.TIF")) # mean annual temperature
ATS <- stack(raster("data_raw/wc2.1_30s_bio_4.TIF")) # temperature seasonality
MAP <- stack(raster("data_raw/wc_MAP_30s.TIF")) # mean annual precipitation
APS <- stack(raster("data_raw/wc2.1_30s_bio_15.TIF")) # precipitation seasonality

# extract the data at the given coordinates
coord <- sites |> dplyr::select(code, latitude, longitude)

pointCoordinates <- sites %>% 
  dplyr::select(code, longitude, latitude) #%>%
 # column_to_rownames("code")

# Transform as coordinates
coordinates(pointCoordinates)= ~ longitude + latitude

MAP_BC<- extract(MAP, pointCoordinates)
MAT_BC <- extract(MAT, pointCoordinates)
ATS_BC <- extract(ATS, pointCoordinates)
APS_BC <- extract(APS, pointCoordinates)

sites <- cbind(sites, MAT_BC = MAT_BC[,1], MAP_BC = MAP_BC[,1],
               ATS_BC = ATS_BC[,1], APS_BC = APS_BC[,1])

#  make them spatial and reproject onto EPSG:3035
sites <- st_as_sf(sites, coords=c('longitude','latitude'), crs=4326)
sites <- st_transform(sites, 3035)

# identify the ecoregion for each site 
sites <- st_join(sites, subset(ecoregions, select=ECO_NAME))

# The GPS for one Estonian site is in the sea, so correct that
sites$ECO_NAME[sites$code == '59_2'] <- "Sarmatic mixed forests"
sites <- sites |> 
  rename(ecoregion = ECO_NAME) |> 
  tibble() |>
  dplyr::select(-geometry) |>
  left_join(coord)

### Load the OTU table and remove OTU not present in any of the remaining plots
OTU_raw <- read.csv("data_raw/OTU.csv") 
  # rename site 50_13 which is wrong into 50_19 in accordance to site data
OTU_raw <- OTU_raw |> mutate(X = ifelse(X == "50_13", "50_19", X))

OTU <- OTU_raw %>%
  select_if(!str_detect(colnames(.), 'NA_SH')) %>%
  column_to_rownames("X")
  
# TRAITS
traits <- read.csv("data_raw/traits.csv",
                   header = TRUE, stringsAsFactors = FALSE)|>
  rename(fruitbody = 'epi..hypogeous', type = simple_type_liberal)


# TAXONOMY
taxa <- read.csv(file = "data_raw/taxonomy.csv", 
                 header = TRUE, stringsAsFactors = FALSE)

# Fix some transcription errors 
taxa$OTU <- sub('Amphimena', 'Amphinema', taxa$OTU)
taxa$OTU <- sub('Elaphomyces_muricatusp', 'Elaphomyces_muricatus__sp', taxa$OTU)

# 3 species with hyphens in binomials that have been turned to dots
taxa$OTU <- sub('(_[a-z]+)\\.([a-z]+_)', '\\1-\\2', taxa$OTU)

# now combine traits and taxonomy to give a single table to filter OTUs
taxa <- merge(taxa, traits, by='OTU', all.x=TRUE)
any(is.na(taxa$fruitbody))

# There are different taxon naming strategies for the taxa/trait data 
# and the OTU table, so reconcile using the SH and denovo part of the tags.
# copy OTU data to modify
otu_regex <- regexpr('denovo_[0-9]+$|SH[0-9]+.07FU$', colnames(OTU))
otu_SH <- regmatches(colnames(OTU), otu_regex)
any(is.na(otu_SH))
length(otu_SH)
colnames(OTU) <- otu_SH

taxa_regex <- regexpr('^denovo_[0-9]+|^SH[0-9]+.07FU$', taxa$SH_code)
taxa$SH <- regmatches(taxa$SH_code, taxa_regex)

otu_to_taxa <- match(colnames(OTU), taxa$SH)
any(is.na(otu_to_taxa))
taxa <- taxa[otu_to_taxa,] 
colnames(OTU) %in% taxa$SH

##### Host specificity ####
# Matrix of occurence
OTU_bin <- OTU
OTU_bin[OTU_bin > 0] <- 1

# matrix OTU per host tree 
OTU_host <- sites |>
  dplyr::select(code, tree_species) |> 
  left_join(OTU_bin |> rownames_to_column("code")) |>
  group_by(tree_species) |>
  summarise_if(is.numeric, sum) |>
  column_to_rownames("tree_species")

OTU_host[OTU_host > 0] <- 1

# Add column of specificity 
OTU_specificity <- as_tibble(cbind(SH = names(OTU_host), t(OTU_host))) |>
  mutate(BEECH = as.numeric(BEECH), 
         OAK = as.numeric(OAK),
         PINE = as.numeric(PINE),
         SPRUCE = as.numeric(SPRUCE)) |>
  rowwise() |>
  mutate(specificity = case_when(
    sum(PINE, SPRUCE) > 0 & sum(BEECH, OAK) > 0  ~ 'generalist',
    sum(PINE, SPRUCE) == 2 & sum(BEECH, OAK) == 0 ~ 'conifer',
    sum(PINE, SPRUCE) == 0 & sum(BEECH, OAK) == 2 ~ "broadleaf",
    PINE == 1 & sum(BEECH, OAK, SPRUCE) == 0 ~ "pine",
    SPRUCE == 1 & sum(BEECH, OAK, PINE) == 0 ~ "spruce",
    BEECH == 1 & sum(SPRUCE, OAK, PINE) == 0 ~ "beech",
    OAK == 1 & sum(SPRUCE, BEECH, PINE) == 0 ~ "oak"
  )) |>
  dplyr::select(SH, specificity)

trait_phylo <- taxa |> 
  left_join(OTU_specificity) |> 
  dplyr::select(SH, phylum, class, order, family, genus, fruitbody, type, specificity, OTU)
colnames(OTU) %in% trait_phylo$SH

## Data selection ----
# Keep ecoregions with more than 3 plots and keep species that are present 
# in at least 1 plot
n_ecoregion <- sites |>  
  count(ecoregion)|>
  filter(n >= 3) # minimum of sites per ecoregion 

## Filter to keep ecoregions with at least 3 plots
sites <- sites |>
  filter(ecoregion %in% n_ecoregion$ecoregion)

### remove OTU not present in any of the remaining plots
OTU <- OTU %>%
  #rownames_to_column("code") %>%
  filter(rownames(OTU) %in% sites$code) %>%
  #column_to_rownames("code")%>%
  # remove OTU with no individuals in the selected sites
  select_if(colSums(.) != 0) %>% #129 sites and 1350 OTUs
  rownames_to_column("code")

# still some species in trait dataset that are not present in the OTU table. Filter.
sp_list <- colnames(OTU)[-1] # remove code name
trait_phylo <- trait_phylo |> 
  filter(SH %in% sp_list)

write_csv(trait_phylo, "data/trait_phylo.csv")
write_csv(OTU, "data/OTU.csv")
write_csv(sites, "data/sites.csv")

#### Prepare a RDS file with all the data 
trait <- read.csv("data/trait_phylo.csv")
OTU <- read.csv("data/OTU.csv")
sites <- read.csv("data/sites.csv")

## change GPS coordinates for the site in accordance with ICP guidelines 
sites <- sites |> 
  mutate(
    latitude  = round(floor(latitude) + (latitude - floor(latitude))/100*60, 2), 
    longitude = round(floor(longitude) + (longitude - floor(longitude))/100*60, 2))

data <- list(OTU, sites, trait)
saveRDS(data, "data/data.rds")
