#

# Script to calculate distance of reef check sites to closest river mouth ----


# Libraries
library(sp)
library(sf)
library(rgdal)
library(raster) 
library(terra)
library(RColorBrewer)
#library(geobr)
library(ggspatial)
library(ggrepel)
library(forcats)
library(beepr)
library(here)


# directories ----
w.dir <- here()
sp.dir <- here('spatial')
d.dir <- here('data')


# Load river mouth points ----
dir(sp.dir)
rm <- st_read(paste(sp.dir, "flowpoints_UTMzone10_NHD_and_CaliStreams_joined.shp", sep ='/'))
plot(rm)
glimpse(rm)

# project --

rm2 <- st_transform(rm, crs = 'epsg:4326')

# Load reef check NC locations ----
dir(d.dir)

rc.sites1 <- read.csv(paste(d.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()

site.names <- rc.sites1$site_name

rc.sites <- st_as_sf(rc.sites1, coords = c('longitude', 'latitude'), crs = 4326)
plot(rc.sites)


# get nearest point ----

near <- st_distance(rc.sites, rm2)
near
class(near)

near2 <- as.data.frame(near)
class(near2)
head(near2)
glimpse(near2)
row.names(near2)
row.names(near2) <- site.names
row.names(near2)
names(near2)

# find closest distance to river mouth ----
min.dist <- apply(near2, 1, FUN = min)
min.dist
class(min.dist)
min.dist2 <- as.data.frame(min.dist) %>% glimpse()
min.dist2$site_name <- site.names
glimpse(min.dist2)

min.dist3 <- min.dist2 %>%
  relocate(site_name, .before = min.dist) %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  rename(min_dist_riv_mth = min.dist) %>%
  glimpse()


# add to data ----
dir(d.dir)

data.nc <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  glimpse()

length(levels(data.nc$site_name))  # 25

# join --

data.nc2 <- data.nc %>%
  left_join(min.dist3, by = 'site_name') %>%
  glimpse()


# save ----
write.csv(data.nc2, paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp_rivmth.csv", sep ='/'), row.names = F)
  
  
  
  
  
  
  
  
  