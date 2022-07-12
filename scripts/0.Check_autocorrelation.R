###
### Script created by : Anita Giraldo on 11 April 2022
### Script last updated by : Anita Giraldo on 11 April 2022

## This script checks for spatial autocorrelation in LANDSAT based kelp data --

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

# Clear environment ----
rm(list=ls())

### Set directories ----
#w.dir<- dirname(rstudioapi::getActiveDocumentContext()$path)
m.dir <- here()
ls.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Kelp_Landsat/data" # library of extracted ls data
s.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Spatial_data/shapefiles"


## Load kelp data ----

df <- read.csv(paste(ls.dir, "NC_Landsat_kelp_area_1984_2021.csv", sep ='/')) %>%
  mutate_at(vars(year_quarter, year, quarter), list(as.factor)) %>%
  glimpse()

# remove Nas --
df2 <- na.omit(df)

# subsample sp --
df3 <- sample_n(df2, 10000)

# convert simple data frame into a spatial data frame object
dfsp <- df3
coordinates(dfsp) <-  ~lon+lat
proj4string(dfsp) <- "+proj=longlat +datum=WGS84 +no_defs"

# Transform sp --
crs1 <- CRS('+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +ellps=GRS80 +units=m +no_defs')
dfsp2 <- spTransform(dfsp, crs1)

# compute variogram --
# https://gsp.humboldt.edu/olm/R/04_01_Variograms.html
TheVariogram <- variogram(area~1, data=dfsp2)
plot(TheVariogram)
TheVariogram


TheVariogramModel <- vgm(psill=0.15, model="Gau", nugget=0.0001, range=5)
plot(TheVariogram, model=TheVariogramModel) 

FittedModel <- fit.variogram(TheVariogram, model=TheVariogramModel)    
plot(TheVariogram, model=FittedModel)


###

## Another way ----
# https://stats.oarc.ucla.edu/r/faq/how-do-i-generate-a-variogram-for-spatial-data-in-r/
# https://zia207.github.io/geospatial-r-github.io/semivariogram-modeling.html

head(df3)

dist1 <- dist(df3[,1:2])
summary(dist1)

v1 <- variog(dfsp, dfsp@coords, dfsp$area)

breaks = seq(0, 1.5, l = 11)
v1

plot(v1, type = "b", main = "Variogram: Av8top")

###

## Another way ----

v.cloud <- variogram(area ~ 1, data = dfsp, cloud=F)
head(v.cloud)
plot(v.cloud, main = "Variogram - default", xlab = "Separation distance (m)")

v.cut <- variogram(area ~ 1, dfsp, cutoff = 400, width = 50)
plot(v.cut, main = "Variogram with cutoff and fixed width", xlab = "Separation distance (m)", xlim = c(0,400))

## fit model 
# show.vgms()
m.exp <- vgm(200, "Spl", 290, 4000, kappa = 0.5)
m.exp <- vgm(200, "Wav", 290, 4000, add.to = vgm(1000, "Exp", 60, nugget = 2.5))
m.exp <- vgm(250, "Wav", 270, 4000, add.to = vgm(200, "Exp", 60, nugget = 10000))
m.exp1 <- vgm(9000, "Wav", range = 270, nugget = 4000)
m.exp2 <- vgm(9000, "Wav", range = 270, nugget = 3000)
m.exp3 <- vgm(10000, "Wav", range = 270, nugget = 4000)
# least square fit
m.exp.f <- fit.variogram(v.cut, m.exp2)
m.exp.f

m.exp3 <- vgm(4952.881, "Wav", range = 0, nugget = 4000)
m.exp3 <- vgm(4952.881, "Nug", range = 0, add.to = vgm(5747.225, "Wav", range = 271.8793))

m.exp.f <- fit.variogram(v.cut,  m.exp3)
m.exp.f


plot(v.cut, pl=F, model=m.exp.f,col="black", cex=1, lwd=0.5,lty=1,pch=20, #xlim = 300,
     main="Variogram models - Wave",xlab="Distance (m)",ylab="Semivariance")

# goodness of fit --
attributes(m.exp.f)$SSErr
