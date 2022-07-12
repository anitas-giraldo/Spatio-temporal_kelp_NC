
## SCRIPT TO EXTEND WAVE DATA -----


library(terra)
library(raster)
library(sf)
library(sp)
library(rgeos)
library(rgdal)
library(class)
library(Rfast)
library(dplyr)
library(Biobase) # matchpt function


###

## load depth data ----

depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Depth"

d <- raster(paste(depth.dir, "depth_allCA_wInterp_wCIN_1km_30depth_latlon.tif", sep ='/'))
plot(d)
d


## load wave data ----

w.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Wave_data"
w1.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Waves_Mainland"

## Max Wave height ----

wh <- stack(paste(w1.dir, "WaveYear_Mean.tif", sep ='/'))
wh

## make depth template ----
d
d[d < 0] <- 1
plot(d)

d.mat <- as.matrix(d)
d.df <- raster::as.data.frame(d, xy = T)

wave.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Wave_data/Wave_extended2"


# LOOP ----
variable <- "wy_mean"
w.years <- paste(2004:2021)

for(i in 1:18) {

### get waves ----

w <- wh[[i]]
plot(w)

w[w == 0] <- NA
plot(w)


w.mat <- as.matrix(w)
w.df <- raster::as.data.frame(w, xy = T)


#### MATCH points ----

## Remove NAs from data
d.df2 <- na.omit(d.df)
w.df2 <- na.omit(w.df)

# Check how many wave data cells ----
w.cells <- length(w.df2$x)

w.df2$index <- paste(1:w.cells)
head(w.df2)


# transform df to matrices ----
d.v <- as.matrix(d.df2)
w.v <- as.matrix(w.df2)


# match points ----
matched <- matchpt(d.v[,1:2], w.v[,1:2])
length(matched$index)

# add coordinates to depth points ----
d.df3 <- cbind(d.df2[,1:2], matched)
head(d.df3)
glimpse(d.df3)


## Join to depth data ----
head(w.df2)
w.df2$index <- as.integer(w.df2$index)

df.join <- d.df3 %>%
  left_join(w.df2, by = "index") %>%
  glimpse()


spg <- df.join
head(spg)
spg <- spg[,c(1:2,7)]

coordinates(spg) <- ~ x.x + y.x
gridded(spg) <- TRUE
rasterDF <- raster(spg)
rasterDF
plot(rasterDF)
proj4string(rasterDF) <- proj4string(w)

raster_name <- paste(paste(variable, w.years[i], sep = '_'), "tif", sep = '.')

writeRaster(rasterDF, paste(wave.dir, variable, raster_name, sep ='/'), overwrite = T)

}







###

### 

###
 
## OPTION 1 tried, using focal analysis ----




fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),2) )
  } else {
    return( round(x[i],2) )
  }
}  

# LOOP ----

w.years <- paste(2004:2021)

for (i in 1:18){
  
  r <- wh_max[[i]]
  
  r2 <- focal(r, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  
  writeRaster(r2, paste(w.dir, "Wave_extended", paste(w.years[i], "extended_nc_wh_max.tif", sep='_'), sep ='/'))
  
}


## Stack wh_max ----

# load raster data --
w2.dir <- paste(w.dir, "Wave_extended", sep ='/')

n.files <- dir(w2.dir)
# list files in source --
n.files <- list.files(w2.dir, pattern = 'wh_max.tif', full.names = TRUE)
n.files
length(n.files)
# list names to load onto the Environment --
names.list <- list.files(w2.dir, pattern = 'wh_max.tif')
names.list <- str_replace(names.list, "wh_max.tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
w.stack <- c()

# use do call to create a raster otherwise it creates a list
w.stack <- do.call("c", tfiles)

# write stack
writeRaster(w.stack, paste(w.dir, "Wave_extended", "extended_nc_wh_max_stack.tif", se='/'))

###

###

###

### WAVE MEANYEAR  ----

w.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Wave_data"

## Max Wave height --

wh_wy <- stack(paste(re.dir, "Waves_mainland",  "WaveYear_Mean.tif", sep ='/'))
wh_wy
plot(wh_wy[[1]])


fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),2) )
  } else {
    return( round(x[i],2) )
  }
}  

# LOOP ----

w.years <- paste(2004:2021)

for (i in 1:18){
  
  r <- wh_wy[[i]]
  
  r2 <- focal(r, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  r2 <- focal(r2, w = matrix(1,3,3), fun = fill.na, 
              pad = TRUE, na.rm = FALSE )
  
  
  writeRaster(r2, paste(w.dir, "Wave_extended", paste(w.years[i], "extended_nc_mean_wy.tif", sep='_'), sep ='/'))
  
}


## Stack wh_max ----

# load raster data --
w2.dir <- paste(w.dir, "Wave_extended", sep ='/')

n.files <- dir(w2.dir)
# list files in source --
n.files <- list.files(w2.dir, pattern = 'mean_wy.tif', full.names = TRUE)
n.files
length(n.files)
# list names to load onto the Environment --
names.list <- list.files(w2.dir, pattern = 'mean_wy.tif')
names.list <- str_replace(names.list, "mean_wy.tif", "")
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
w.stack <- c()

# use do call to create a raster otherwise it creates a list
w.stack <- do.call("c", tfiles)

# write stack
writeRaster(w.stack, paste(w.dir, "Wave_extended", "extended_nc_mean_wy_stack.tif", se='/'))
