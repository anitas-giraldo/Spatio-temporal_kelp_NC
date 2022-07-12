###
### Script created by : Anita Giraldo on 11 April 2022
### Script last updated by : Anita Giraldo on 11 April 2022

## This script removes data points that are at a certain distance to address autocorrelation --

# libraries ----
library(spdep)
library(sp)
library(sf)
library(rgdal)
library(rgeos)
library(gstat)
library(here)
library(dplyr)
library(geoR)
library(spThin) # only for occurrence points
library(maptools)
library(spatstat.geom)
library(spatstat.data)
library(spatstat)
library(caret)

# Clear environment ----
rm(list=ls())

### Set directories ----
#w.dir<- dirname(rstudioapi::getActiveDocumentContext()$path)
m.dir <- here()
ls.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Kelp_Landsat/data" # library of extracted ls data
s.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Spatial_data/shapefiles"
d.dir <- here::here('data')


## Load kelp data ----

df <- read.csv(paste(ls.dir, "NC_Landsat_kelp_area_1984_2021.csv", sep ='/')) %>%
  mutate(latlon = paste(lat, lon, sep = '_')) %>%
  mutate_at(vars(year_quarter, year, quarter, latlon), list(as.factor)) %>%
  glimpse() # Rows: 5,748,336

# get maximum area of the year ----
df1 <- df %>%
  group_by(year, latlon) %>%
  summarize(max.area = max(area, na.rm = TRUE)) %>%
  mutate(latitude = sub("_.*", "", latlon)) %>% # Extract characters before pattern _
  mutate(longitude = sub(".*_", "", latlon)) %>% # Extract characters after pattern _
  na_if(-Inf) %>%
  ungroup() %>%
  glimpse() # Rows: 1,437,085

class(df1)

# remove Nas --
df2 <- na.omit(df1) %>%
  mutate_at(vars(latitude, longitude), list(as.numeric)) %>%
  glimpse() # Rows: 1,435,851

any(is.na(df2))

# convert simple data frame into a spatial data frame object
# dfsp <- df2
# coordinates(dfsp) <-  ~lon+lat
# proj4string(dfsp) <- "+proj=longlat +datum=WGS84 +no_defs"

# # Transform sp --
# crs1 <- CRS('+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +ellps=GRS80 +units=m +no_defs')
# dfsp2 <- spTransform(dfsp, crs1)

###

### now I need to remove records closer than 200 m ----
# run the function 'filterByProximity' for every year 

## Filter by proximity function ----
# https://stackoverflow.com/questions/22051141/spatial-filtering-by-proximity-in-r
filterByProximity <- function(xy, dist, mapUnits = F) {
  #xy can be either a SpatialPoints or SPDF object, or a matrix
  #dist is in km if mapUnits=F, in mapUnits otherwise 
  if (!mapUnits) {
    d <- spDists(xy,longlat=T)
  }
  if (mapUnits) {
    d <- spDists(xy,longlat=F)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    return(xy[-discard[,1],])
  }
  if (nrow(closePts) == 0) {
    return(xy)
  }
}

##

## other info for loop ----


glimpse(df2)

length(levels(df2$year)) # 38
years.list <- paste(1999:2021)
crs1 <- CRS('+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +ellps=GRS80 +units=m +no_defs')
crs2 <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# blank df to add years --
all.years.df <- data.frame()

## loop ----

for(i in 1:length(years.list)){
  
  df.year.x <- df2 %>% dplyr::filter(year == years.list[i]) #%>%glimpse()
  
  # first separate zeros and others 
  df.zeros <- df.year.x %>%
    dplyr::filter(max.area == 0) %>% 
    glimpse() # Rows: 35,635
  
  df.no.zero <- df.year.x %>%
    dplyr::filter(max.area > 0) %>% 
    glimpse() # Rows: 2,138
  
  # second remove points closer than 200 m
  dfsp.x <- SpatialPointsDataFrame(coords      = df.no.zero[,c(5,4)],
                                   data        = df.no.zero, 
                                   proj4string = crs2)
  dfsp.x2 <- spTransform(dfsp.x, crs1)
  plot(dfsp.x2$max.area)
  
  test <- filterByProximity(dfsp.x, 0.2, mapUnits = F)
  #plot(test$max.area)
  test.df <- as.data.frame(test)
  
  # third subsample not zeros
  size.1 <- 175/length(test.df$max.area)
  year.sub1 <- createDataPartition(test.df$max.area, p = size.1, 
                                    list = FALSE, 
                                    times = 1)
  year.sub.df1 <- test.df[year.sub1,] #%>%glimpse()
  plot(year.sub.df1$max.area)
  
  # fourth subsample zeros
  #size.2 <- 50/length(df.zeros$max.area)
  year.sub2 <- sample_n(df.zeros, 30)
  
  dfsp.xz <- SpatialPointsDataFrame(coords      = year.sub2[,c(5,4)],
                                   data        = year.sub2, 
                                   proj4string = crs2)
  dfsp.xz2 <- spTransform(dfsp.xz, crs1)
  #plot(dfsp.xz)
  testz <- filterByProximity2(dfsp.xz2, 200)
  #plot(testz)
  test.dfz <- as.data.frame(testz)
  
  # Join zeros and non zeros ---
  glimpse(test.dfz)
  glimpse(year.sub.df1)
  
  df.all.year <- rbind(year.sub.df1, test.dfz)
  
  # add to main df --
  all.years.df <- rbind(df.all.year, all.years.df)
  
  # system.time(test3 <- filterByProximity2(dfsp.x2, dist=100))
  # plot(test)
  # plot(test$max.area)
  # writeOGR(test3, d.dir, "test3", driver = "ESRI Shapefile")
  
}

beep()

## save ----
write.csv(all.years.df, paste(d.dir, "NC_subsampled_data_Landsat_1984-2021.csv", sep ='/'), row.names = F)

beep()


glimpse(all.years.df)

###


## Another proximity function ----
# https://stackoverflow.com/questions/22051141/spatial-filtering-by-proximity-in-r
filterByProximity2 <- function(occ, dist) {
  pts <- as.ppp.SpatialPoints(occ)
  d <- spatstat.geom::nndist(pts)
  z <- which(d > dist)
  return(occ[z,])
}

# x <- filterByProximity2(dfsp, 0.2) ## example
