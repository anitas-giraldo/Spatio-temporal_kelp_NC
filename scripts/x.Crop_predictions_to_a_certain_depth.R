###
### Script based on one created by : Anita Giraldo on 26 April 2022
### Script last updated by : Anita Giraldo on 26 April 2022

## This script crops the predictions of kelp and urchin densitites to a maximum depth of 30 m --


# libraries ----
library(sp)
library(sf)
library(raster)
library(terra)
library(rdgal)
library(rgeos)
library(here)
library(stringr)


# directories ----

m.dir <- here()
ko.dir <- here("outputs_nc_rcca/gam_V3")
k.preds.dir <- paste(ko.dir, "sp_predictions", sep ='/')
uo.dir <- here("outputs_nc_rcca/gam_urchins2")
u.preds.dir <- paste(uo.dir, "sp_predictions")
hd.dir <- "D:/"
sp.dir <- here("spatial")


# 1. Get bathymetry ----

# load bathy at 30m resolution --

bathy <- rast(paste(hd.dir, "CA_agg_bathy_30m_wInterp", "depth_mean_nc.all_wInterp_30m.tif", sep ='/'))
plot(bathy)

# crop to 30m of depth --
bathy30 <- classify(bathy, cbind(-Inf, -31, NA))
plot(bathy30)

# save --
writeRaster(bathy30, paste(sp.dir, "NC_bathy_wInterp_30mdepth_30mres.tif", sep ='/'))

# aggregate to 300m --
fact <- 300/30

bathy30agg <- aggregate(bathy30, fact = fact, fun = mean, na.rm = T)
plot(bathy30agg)

# save --
writeRaster(bathy30agg, paste(sp.dir, "NC_bathy_wInterp_30mdepth_300mres.tif", sep ='/'))


###


# 2. Get kelp predictions ----
dir(k.preds.dir)

# load raster data --
n.files <- dir(k.preds.dir)
# list files in source --
n.files <- list.files(k.preds.dir, pattern = '.tif', full.names = TRUE)
n.files
length(n.files)
# list names to load onto the Environment --
names.list <- list.files(k.preds.dir, pattern = '.tif')
names.list <- str_replace(names.list, ".tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)
plot(preds.stack)


###


# 3. Mask kelp predictions ----

# mask by bathy to a max of 30 m of depth ----
preds30 <- mask(preds.stack, bathy30agg)

# check it worked --
par(mfrow=c(1,2))
plot(preds.stack[[26]])
plot(preds30[[26]])
names(preds30)

# name the layers ----
names.nereo <- paste("Nereo_density", 1998:2021, sep = '_')
names.stack <- c(names.nereo, "mean_nereo" ,    "sd_nereo"   ,    "s_nereo"  , "s_index_nereo")

names(preds30) <- names.stack

## save croped outcome ----
writeRaster(preds30, paste(sp.dir, "Kelp_predictions_gam_V3_30mdepth.tif", sep ='/'), overwrite = T)


### Aggregate stability ----

# 1km ----

s <- preds30[[28]]
s

fact <- 1000/300
s1kmean <- aggregate(s, fact=fact, fun = mean, na.rm = T)
plot(s1kmean)

s1kmax <- aggregate(s, fact=fact, fun = max, na.rm = T)
plot(s1kmax)

s1kmin <- aggregate(s, fact=fact, fun = min, na.rm = T)
plot(s1kmin)


# 2km ----

fact <- 2000/300
s2kmean <- aggregate(s, fact=fact, fun = mean, na.rm = T)
plot(s2kmean)

# 5km ----

fact <- 5000/300
s5kmean <- aggregate(s, fact=fact, fun = mean, na.rm = T)
plot(s5kmean)


### Make 1km spatial points ----

s1kdf <- as.data.frame(s1kmean, xy = T)
head(s1kdf)
any(is.na(s1kdf))

write.csv(s1kdf, paste(sp.dir, "Nereo_stability_index_points_1km.csv", sep ='/'), row.names = F)
