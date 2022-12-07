
## Script by Anita Giraldo, 2 September2022
## Last modified by Anita Giraldo, 5 Dec 2022

##  This script extends the pixels of environmental predictors to encompass all of the
## desired predicted area


# load libraries ----
library(terra)
library(raster)
library(sf)
library(sp)
library(rgeos)
library(rgdal)
library(class)
library(Rfast)
library(dplyr)
library(here)
library(stringr)
library(Biobase) # matchpt function

###

# * Days 10 NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# nitrate dir
nit.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"

## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(nit.dir, "Days_10N.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)

# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Days_10N"

# create directory for outputs --
#dir.create(paste(o2.dir, variable, sep ='/'))

# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))

###

###

###



# * MAX NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# nitrate dir
nit.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"

## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(nit.dir, "Max_Monthly_Nitrate.tif", sep = '/'))
plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)

# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Max_Monthly_Nitrate"

# create directory for outputs --
#dir.create(paste(o2.dir, variable, sep ='/'))

# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))

###

###

###


# * MEAN MONTHLY TEMP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# temp dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Temperature"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Temperature", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Mean_Monthly_Temp"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Mean_Monthly_Temp.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MEAN MONTHLY UPWELLING TEMP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# temp dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Temperature"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Temperature", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Mean_Monthly_Upwelling_Temp"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Mean_Monthly_Upwelling_Temp.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))

###


###


###


# * DAYS 16 C ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# temp dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Temperature"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Temperature", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Days_16C"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Days_16C.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))

###


###


###


# * MIN MONTHLY TEMP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# temp dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Temperature"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Temperature", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_Temp"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_Temp.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###


# * MIN MONTHLY ANOMALY TEMP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# temp dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Temperature"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Temperature", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_Anomaly_Temp"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_Anomaly_Temp.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###


# * MAX MONTHLY ANOMALY SUMMER NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Max_Monthly_Anomaly_Summer_Nitrate"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Max_Monthly_Anomaly_Summer_Nitrate.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###


# * MIN MONTHLY NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_Nitrate"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_Nitrate.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MIN MONTHLY ANOMALY NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_Anomaly_Nitrate"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_Anomaly_Nitrate.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * WH MEAN ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
#env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Waves"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables/Waves"


### THIS WAS ALREADY EXTENDED, BUT I SAVED THE STACKED RASTER --

# set output dir ----
o2.dir <- paste(o.dir, "archive/wh_mean", sep = '/')
o2.dir

# set output dir ----
o3.dir <- paste(o.dir, "wh_mean", sep = '/')
o3.dir

# LOOP ----

# variable to extend --
variable <- "wh_mean"

# create directory for outputs --
#dir.create(paste(o2.dir, variable, sep ='/'))

## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(o2.dir, "wh_mean_extended_stack.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
# 
# d2 <- aggregate(d, fact = 5, fun = 'mean', na.rm =T)
# plot(d2)
# d2

e <- raster::resample(env.rast, d)
plot(e[[1]])

# d3 <- raster::resample(d, e)
# plot(d3)
# d3

d3 <- d

## make depth template ----
d3
d3[d3 < 0] <- 1
plot(d3)


d.mat <- as.matrix(d3)
d.df <- raster::as.data.frame(d3, xy = T)



# get years --
rast.years <- paste(2004:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- e[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MIN MONTHLY ANOMALY SUMMER NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_Anomaly_Summer_Nitrate"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_Anomaly_Summer_Nitrate.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----

d2 <- aggregate(d, fact = 5, fun = 'mean', na.rm =T)
plot(d2)
d2

e <- raster::resample(env.rast, d2)
plot(e[[1]])

# d3 <- raster::resample(d, e)
# plot(d3)
# d3

d3 <- d2

## make depth template ----
d3
d3[d3 < 0] <- 1
plot(d3)


d.mat <- as.matrix(d3)
d.df <- raster::as.data.frame(d3, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- e[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))

###


###


###

# * MAX MONTHLY NPP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/NPP"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "NPP", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Max_Monthly_NPP"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Max_Monthly_NPP.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----

d2 <- aggregate(d, fact = 5, fun = 'mean', na.rm =T)
plot(d2)
d2

e <- raster::resample(env.rast, d2)
plot(e[[1]])

# d3 <- raster::resample(d, e)
# plot(d3)
# d3

d3 <- d2

## make depth template ----
d3
d3[d3 < 0] <- 1
plot(d3)


d.mat <- as.matrix(d3)
d.df <- raster::as.data.frame(d3, xy = T)



# get years --
rast.years <- paste(2003:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- e[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MIN MONTHLY NPP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/NPP"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "NPP", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_NPP"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_NPP.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(2003:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###


# * MAX MONTHLY NPP UPWELLING ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/NPP"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "NPP", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Max_Monthly_NPP_Upwelling"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Max_Monthly_NPP_Upwelling.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(2003:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MAX ANOM UPWELLING NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Max_Monthly_Anomaly_Upwelling_Nitrate"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Max_Monthly_Anomaly_Upwelling_Nitrate.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MIN NPP UPWELLING----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/NPP"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "NPP", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Min_Monthly_NPP_Upwelling"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Min_Monthly_NPP_Upwelling.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(2003:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MAX MONTHLY ANOMALY UPWELLING TEMP ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# temp dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Temperature"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Temperature", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Max_Monthly_Anomaly_Upwelling_Temp"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Max_Monthly_Anomaly_Upwelling_Temp.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----
d <- raster::resample(d, env.rast)


## make depth template ----
d
d[d < 0] <- 1
plot(d)


d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)



# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- env.rast[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))



###


###


###

# * MEAN NITRATE ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Nitrate"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "Nitrate", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "Mean_Monthly_Nitrate"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load nitrate data ----
env.rast <- stack(paste(env.dir, "Mean_Monthly_Nitrate.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])


## resample depth to env variable ----

d2 <- aggregate(d, fact = 5, fun = 'mean', na.rm =T)
plot(d2)
d2

e <- raster::resample(env.rast, d2)
plot(e[[1]])

# d3 <- raster::resample(d, e)
# plot(d3)
# d3

d3 <- d2

## make depth template ----
d3
d3[d3 < 0] <- 1
plot(d3)


d.mat <- as.matrix(d3)
d.df <- raster::as.data.frame(d3, xy = T)




# get years --
rast.years <- paste(1998:2021)
length(rast.years)


for(i in 1:length(rast.years)) {
  
  ### get waves ----
  
  rast.year <- e[[i]]
  plot(rast.year)
  
  # rast.year[rast.year == 0] <- NA
  # plot(rast.year)
  
  
  rast.mat <- as.matrix(rast.year)
  rast.df <- raster::as.data.frame(rast.year, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  proj4string(rasterDF) <- proj4string(rast.year)
  
  raster_name <- paste(paste(variable, rast.years[i], sep = '_'), "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  
}


## Stack the extended variable and save ----

# load raster data --
o3.dir <- paste(o2.dir, variable, sep ='/')

x.files <- dir(o3.dir)
x.files
# list files in source --
x.files <- list.files(o3.dir, pattern = '.tif$', full.names = TRUE)
x.files
length(x.files)
# list names to load onto the Environment --
names.list <- list.files(o3.dir, pattern = '.tif$')
names.list <- str_replace(names.list, ".tif$", "")
length(names.list)
names.list

# load csv files as a list --
rfiles <- lapply(x.files, rast) # this is a list
rfiles[[1]]
plot(rfiles[[1]])

# stack them ---
r.stack <- c()

# use do call to create a raster otherwise it creates a list
r.stack <- do.call("c", rfiles)

# write stack
writeRaster(r.stack, paste(o3.dir, paste(variable, "extended_stack.tif", sep='_'), sep='/'))


###


###


###

# * MEAN VRM ----

# Clear environment ----
rm(list=ls())

# directories ----
m.dir <- here()
# depth dir
depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"
# env dir
env.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"
# output dir
o.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"




# set output dir ----
o2.dir <- paste(o.dir, "VRM", sep = '/')
o2.dir

# LOOP ----

# variable to extend --
variable <- "mean_vrm"

# create directory for outputs --
dir.create(paste(o2.dir, variable, sep ='/'))


## load depth data ----

d <- rast(paste(depth.dir, "depth_allCA_wInterp_wCIN_300res_30m_latlon.tif", sep ='/'))
plot(d)
d

# if need to aggregate --
# d2 <- classify(d, cbind(-Inf, -31, NA))
# d2
# 
# d3 <- aggregate(d2, fact = 10, fun = 'mean', na.rm = TRUE)
# plot(d3)
# 
# d3 <- terra::project(d3, "epsg:4326")
# 
# 
# writeRaster(d3, paste(depth.dir, "depth_allCA_wInterp_wCIN_300res_30m_latlon.tif", sep ='/'), overwrite = T)



## load nitrate data ----
env.rast <- rast(paste(env.dir, "vrm_mean_cc_all-mapped_300m_wInterp.tif", sep = '/'))
#plot(env.rast)
env.rast
plot(env.rast[[1]])

env.rast <- terra::project(env.rast, "epsg:4326")
plot(env.rast)


## resample depth to env variable ----

d2 <- crop(d, env.rast)
plot(d2)

e <- terra::resample(env.rast, d2)
plot(e)

d3 <- terra::resample(d2, e)
plot(d3)



stack1 <- c(e, d3)
plot(stack1)

## make depth template ----
d3
d3[d3 < 0] <- 1
plot(d3)


d.mat <- as.matrix(d3)
d.df <- raster::as.data.frame(d3, xy = T)




  
  rast.mat <- as.matrix(e)
  rast.df <- raster::as.data.frame(e, xy = T)
  
  
  #### MATCH points ----
  
  ## Remove NAs from data
  d.df2 <- na.omit(d.df) # depth df
  rast.df2 <- na.omit(rast.df) 
  
  # Check how many wave data cells ----
  rast.cells <- length(rast.df2$x)
  
  rast.df2$index <- paste(1:rast.cells)
  head(rast.df2)
  
  
  # transform df to matrices ----
  d.v <- as.matrix(d.df2)
  rast.v <- as.matrix(rast.df2)
  
  
  # MATCH POINTS ----
  matched <- matchpt(d.v[,1:2], rast.v[,1:2])
  length(matched$index)
  
  # add coordinates to depth points ----
  d.df3 <- cbind(d.df2[,1:2], matched)
  head(d.df3)
  glimpse(d.df3)
  
  
  ## Join to depth data ----
  head(rast.df2)
  rast.df2$index <- as.integer(rast.df2$index)
  
  df.join <- d.df3 %>%
    left_join(rast.df2, by = "index") %>%
    glimpse()
  
  
  spg <- df.join
  head(spg)
  spg <- spg[,c(1:2,7)]
  
  coordinates(spg) <- ~ x.x + y.x
  gridded(spg) <- TRUE
  rasterDF <- raster(spg)
  rasterDF
  plot(rasterDF)
  crs(rasterDF) <- crs(env.rast)
  
  raster_name <- paste(variable, "tif", sep = '.')
  
  writeRaster(rasterDF, paste(o2.dir, variable, raster_name, sep ='/'), overwrite = T)
  writeRaster(d3, paste(o2.dir, "d_test.tif", sep ='/'), overwrite = T)




