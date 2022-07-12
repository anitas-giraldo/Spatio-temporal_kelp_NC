
## Script by Anita Giraldo, 21 June 2022
## Last modified by Anita Giraldo, 21 June 2022

##  This script predicts and projects spatially the result from the models selected in script 2

# libraries ----
library(Amelia) # to map missing data
library(here)
library(tidyverse)
library(dplyr)
library(ggplot2) 
library(caret) ## For model fitting and evaluation
library(RANN)
library(corrplot)
library(rsample)
library(yardstick) ## to evaluate a model
library(visreg)  ## For visualizing regression models
library(plotROC) ## For constructing ROC curves
library(mgcv)    ## For fitting GAM models
library(kernlab) ## Contains an example dataset
library(glmnet)  ## For fitting regularized models
library(pROC)
library(ROCR)
library(car)
library(OddsPlotty) # to plot from the caret package
library(tidymodels) # modern version of caret
library(cutpointr) # Find optimal cutoff point for binary classification
library(bestglm) # Best subset GLM : The function bestglm selects the best subset of inputs for the glm family
library(MASS) # Stepwise regression
library(reshape2)
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
library(gridExtra)
library(FSSgam)
library(viridis)
library(hrbrthemes)
library(gplots)
library(beepr)
library(stringr)



# Clear environment ----
rm(list=ls())



# directories ----
m.dir <- here()
d.dir <- here('outputs_nc_rcca')
k.dir <- paste(d.dir, "gam_V4", sep ='/') # kelp model results
u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"
o.dir <- here("outputs_nc_rcca/gam_V4/gam_5.1")
sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

###
###

# Model V4 5.1 ####

## Load probability of rock layer ----

dir(sub.dir)

rock <- rast(paste(sub.dir, "prob_rock_nc.all_300m_wInterp.tif", sep = '/')) 
plot(rock)

# load depth upt to 30 m ----

depth <- rast(paste(sub.dir, "depth_mean_nc.all_wInterp_300m_30m.tif", sep ='/'))
plot(depth)

rock.mask <- mask(rock, depth, updatevalue = 1000)
plot(rock.mask)
rock.mask[is.na(rock.mask)] <- 0
rock.mask[rock.mask == 1000] <- NA
plot(rock.mask)

par(mfrow=c(1,3),mar=c(2,4,3,1))
plot(rock)
plot(depth)
plot(rock.mask)
dev.off()

# project rock layer --
crs1 <- "epsg:4326"
rock.mask.p <- project(rock.mask, crs1)
plot(rock.mask.p)


# save Masked Rock projected --
# writeRaster(rock.mask.p, paste(sub.dir, "prob_rock_nc_MASKED_LatLong_all_300m_wInterp.tif", sep = '/'))
# writeRaster(rock.mask, paste(sub.dir, "prob_rock_nc_MASKED_UTM_all_300m_wInterp.tif", sep = '/'))

###

###


## Load spatial predictions  ----

# define folder with predictions
preds.dir <- paste(o.dir, "sp_predictions_v3", sep ='/')
preds.dir

# load raster data --

n.files <- dir(preds.dir)
# list files in source --
n.files <- list.files(preds.dir, pattern = '3.tif', full.names = TRUE)
n.files
length(n.files)
n.files<-n.files[-19]
length(n.files)

# list names to load onto the Environment --
names.list <- list.files(preds.dir, pattern = '3.tif')
names.list <- str_replace(names.list, "3.tif", "")
length(names.list)
names.list <- names.list[-19]
length(names.list)

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]] 

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)


### 

###


## Scale predictions to probability of rock ----

## resample --

rock.mask.p2 <- resample(rock.mask.p, preds.stack[[1]])


## scale --

preds.stack.scaled <- preds.stack*rock.mask.p2

par(mfrow=c(1,2),mar=c(2,4,3,1))
plot(preds.stack[[10]])
plot(preds.stack.scaled[[10]])
dev.off()

# save 
#writeRaster(preds.stack.scaled, paste(preds.dir, "Sp_predictions_Nereo_NC_gam_V4_5.1_v3_Standarized.tif", sep ='/'), overwrite = T)


###

###

## Recalculate stability ----

stability.stack <- c()

# calculate mean across years
mean.nereo <- mean(preds.stack.scaled, na.rm = T)
plot(mean.nereo)


# calculate sd across years
sd.nereo <- app(preds.stack.scaled, fun = sd, na.rm = T)
plot(sd.nereo)


# calcualte standarized mean

# get maximum mean
max.mean <- minmax(mean.nereo)[2] 

st.mean <- mean.nereo/max.mean
plot(st.mean)


# calculate stability
s.nereo <- mean.nereo/sd.nereo
plot(s.nereo)

# calculate stability index
s_index.nereo <- (s.nereo*mean.nereo)
plot(s_index.nereo)


# calculate stability index with sandarised mean
s_index2.nereo <- (s.nereo*st.mean)
plot(s_index2.nereo)


# stack 

stability.stack <- c(mean.nereo, sd.nereo, st.mean, s.nereo, s_index.nereo)
names(stability.stack) <- c("mean_nereo", "sd_nereo", "st_mean", "s_nereo", "s_index_nereo")
plot(stability.stack)

# save 
#writeRaster(stability.stack, paste(preds.dir, "Stabilty_calculations_Nereo_NC_gam_V4_5.1_v3_Standarized.tif", sep ='/'), overwrite = T)


###

###

## Save as csv to plot as points ----


## Mean --
mean.points <- as.data.frame(mean.nereo, xy = T)
head(mean.points)

# save
write.csv(mean.points, paste(preds.dir, "mean_nereo_NC_gam_V4_5.1_v3_Standarized.csv", sep ='/'))


## Stability  --
s.points <- as.data.frame(s.nereo, xy = T)
head(s.points)
names(s.points) <-  c("x"  ,  "y"  ,  "stability")

# save
write.csv(s.points, paste(preds.dir, "stability_nereo_NC_gam_V4_5.1_v3_Standarized.csv", sep ='/'))


## Stability  index --
s_index.points <- as.data.frame(s_index.nereo, xy = T)
head(s_index.points)
names(s_index.points) <-  c("x"  ,  "y"  ,  "stability_index")

# save
write.csv(s_index.points, paste(preds.dir, "stability_index_nereo_NC_gam_V4_5.1_v3_Standarized.csv", sep ='/'))
