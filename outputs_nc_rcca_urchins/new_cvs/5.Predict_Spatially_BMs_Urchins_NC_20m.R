##

## Script by Anita Giraldo, 4 May 2022
## Last modified by Anita Giraldo, 6 Dec 2022

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



# clear environment ----
rm(list = ls())


# directories ----
m.dir <- here()
d.dir <- here('data')
o.dir <- here('outputs_nc_rcca_urchins')
#k.dir <- paste(o.dir, "gam_urchins4", sep ='/')
cv.dir <- paste(o.dir, "new_cvs", sep ='/')
#k.dir <- paste(o.dir, "gam_urchins5", sep ='/') # kelp model results
#u.dir <- paste(d.dir, "gam_urchins3", sep ='/') # urchin model results
#rcca.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/North_Coast_w_RCCA/raw_data"
#dd.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/R_Projects/Extract_env_data/nc.rcca.outputs"



# 1. PREPARE DATA ####

### 1.1. Load info on years RCCA ----
years <- read.csv(paste(d.dir, "RCCA_North_Coast_sites.csv", sep ='/')) %>%
  glimpse()


### 1.2. Get sites with preMHW data ----
# 3 or more pre MHW surveys
ncsites <- years %>%
  mutate_at(vars(site_name), list(as.factor)) %>%
  # get only sites with PRE MHW data 
  dplyr::filter(pre.mhw.years > 2) %>%
  droplevels() %>%
  glimpse() # Rows: 10


### 1.3. Load RCCA data ----

df <- read.csv(paste(d.dir, "RCCA_kelp_inverts_NC_depth-zones_wave_clim_temp_nit_subs_orbvel_npp.csv", sep ='/')) %>%
  mutate_at(vars(site_name, month, year, transect, zone), list(as.factor)) %>%
  glimpse() # Rows: 1,154


### 1.4. Get the sites for North Coast model ----
# sites selected in 'ncsites'
df.nc <- df %>%
  dplyr::select(-c(latitude, longitude)) %>%
  right_join(ncsites, by = c('site_name')) %>%
  droplevels() %>% #glimpse()
  dplyr::select(-c(total.years, pre.mhw.years, during.mhw.years, post.mhw.years)) %>%
  relocate(c(latitude, longitude), .after = zone) %>%
  glimpse() # Rows: 708

length(levels(df.nc$site_name)) # 10
levels(df.nc$site_name)
any(is.na(df.nc$Max_Monthly_Anomaly_Temp))


### 1.5. Choose variables and transform needed ----
# as per the FSSgam
names(df.nc)

dat1 <- df.nc %>%
  dplyr::select(
    # Factors 
    latitude, longitude,
    site_name, year, transect, zone,
    # Bio vars
    den_NERLUE , den_MESFRAAD , den_STRPURAD , den_PYCHEL, den_HALRUF,
    # Nitrate vars 
    Days_10N, 
    Min_Monthly_Nitrate, 
    Max_Monthly_Nitrate,
    Mean_Monthly_Nitrate,
    Mean_Monthly_Upwelling_Nitrate,
    Max_Monthly_Anomaly_Nitrate,
    Mean_Monthly_Summer_Nitrate,
    # Temperature vars
    Days_16C ,
    Mean_Monthly_Temp ,
    Mean_Monthly_Summer_Temp,
    MHW_Upwelling_Days  , 
    Min_Monthly_Anomaly_Temp,
    Max_Monthly_Anomaly_Upwelling_Temp,
    Min_Monthly_Temp, 
    Mean_Monthly_Upwelling_Temp,
    # Other vars
    #npp.mean,                             
    #wh.95 ,   wh.max,
    npgo_mean , mei_mean,
    # substrate
    mean_depth, mean_prob_of_rock, mean_vrm, mean_slope,
    # waves
    wh_max, wh_mean, mean_waveyear, wh_95prc,
    # Orb vel
    UBR_Mean, UBR_Max,
    # NPP
    Mean_Monthly_NPP, Max_Monthly_NPP_Upwelling, Mean_Monthly_NPP_Upwelling, Min_Monthly_NPP) %>%
  # Bio transformations
  mutate(log_den_NERLUE = log(den_NERLUE + 1),
         log_den_MESFRAAD = log(den_MESFRAAD + 1),
         log_den_STRPURAD = log(den_STRPURAD + 1),
         log_den_PYCHEL = log(den_PYCHEL + 1),
         log_den_HALRUF = log(den_HALRUF + 1),
         log_mean_vrm = log(mean_vrm + 1)) %>%
  dplyr::select(-c(den_NERLUE,
                   den_MESFRAAD,
                   den_STRPURAD,
                   den_PYCHEL,
                   den_HALRUF,
                   mean_vrm)) %>%
  # Temperature transformations
  mutate(log_Days_16C = log(Days_16C + 1)) %>%
  dplyr::select(-c(Days_16C)) %>%
  # Orb vel transformations
  mutate(log_UBR_Mean = log(UBR_Mean + 1),
         log_UBR_Max = log(UBR_Max + 1)) %>%
  dplyr::select(-c(UBR_Mean,
                   UBR_Max)) %>%
  # NPP transformations
  mutate(log_Mean_Monthly_NPP_Upwelling = log(Mean_Monthly_NPP_Upwelling + 1),
         log_Min_Monthly_NPP = log(Min_Monthly_NPP + 1)) %>%
  dplyr::select(-c(Mean_Monthly_NPP_Upwelling,
                   Min_Monthly_NPP)) %>%
  glimpse() # Rows: 708


#### Drop NAs ----
dat2 <- dat1 %>%
  drop_na() %>%
  glimpse() # Rows: 686


glimpse(dat2)
levels(dat2$year)
names(dat2)
min(dat2$mean_depth)

# 2. LOAD BEST MODELS ####


### 2.1. Load best model csv ----
best_mods <- read.csv(paste(cv.dir, "best_model.csv", sep ='/')) 
head(best_mods)

names(best_mods)
names(best_mods) <- c("bm_id"  , "modname")
names(best_mods)



### 2.2. Set submodel names ----
no_submodels <- nrow(best_mods)
no_submodels


#submodels <- c('3.0', '3.1', '3.2', '3.3')
submodels <- paste(1:no_submodels)

submodel_names <- paste("bm", submodels, sep = '')



# 3.1. Get submodel name --
bm_name <- best_mods$modname[1]

# get submodel variables --
bm_vars <- strsplit(bm_name, split = "+", fixed = T)

# 3.2. Set model formula --

# set dependent variable --
dep <- 'log_den_STRPURAD'

# set predictors --
preds <- bm_vars[[1]]

# get variables with smoothers --

preds2 <- character()

for (j in 1:length(preds)) {
  pred.x <- preds[j]
  pred_form <- paste("s(", pred.x, ", k = 4, bs = 'cr')", sep = '')
  preds2 <- append(preds2, pred_form)
}

preds2


# get random variables --
random_vars <- "s(site_name, zone, bs = 're') + s(year, bs = 're')"


# set model formula --
bm_form <- as.formula(paste(dep, paste(paste(preds2, collapse = " + "), random_vars, sep = "+"), sep=" ~ "))
bm_form



# Fit model bm 3 using all data ----
#bm1 <- gam(bm_form, family = tw(), data = dat2, method = "GCV.Cp")
bm <- gam(bm_form, family = tw(), data = dat2, method = "REML")

bm$sp[2]

# check model ----

bm$aic
summary(bm)
gam.check(bm)
summary(bm)$dev.expl

par(mfrow=c(3,3),mar=c(2,4,3,1))
visreg(bm)
dev.off()





###


###


# PREDICT BEST MODEL V2 - bm3 plus max npp ----

# directory of extended variables ----
re.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Extended_variables"


### Get slope ----

depth.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

dir(depth.dir)

preds

d1 <- rast(paste(depth.dir, "depth_mean_nc_all_wInterp_300m_30m.tif", sep ='/'))
plot(d1)
d1

d <- classify(d1, cbind(-Inf, -20, NA))
plot(d)

depth <- terra::project(d, "epsg:4326")
plot(depth)


n.extent <- ext(depth)



### Get NPP predictors  ----

# Mean_Monthly_NPP_Upwelling ----

npp.up <- rast(paste(re.dir, 'NPP', "Mean_Monthly_NPP_Upwelling", "Mean_Monthly_NPP_Upwelling_extended_stack.tif", sep ='/'))
npp.up
plot(npp.up[[2]])



# crop to NC --
npp.up2 <- crop(npp.up, n.extent)
plot(npp.up2[[1]])

# resample predictors to bathy --
npp.up3 <- resample(npp.up2, depth)
plot(npp.up3[[2]])


# mask predictors to bathy --
npp.up4 <- mask(npp.up3, depth)
plot(npp.up4[[2]])

# log --
npp.up5 <- log(npp.up4 + 1)
plot(npp.up5[[2]])


# Min_Monthly_NPP ----

npp <- rast(paste(re.dir, 'NPP', "Min_Monthly_NPP", "Min_Monthly_NPP_extended_stack.tif", sep ='/'))
npp
plot(npp[[2]])



# crop to NC --
npp2 <- crop(npp, n.extent)
plot(npp2[[1]])

# resample predictors to bathy --
npp3 <- resample(npp2, depth)
plot(npp3[[2]])


# mask predictors to bathy --
npp4 <- mask(npp3, depth)
plot(npp4[[2]])

# log --
npp5 <- log(npp4 + 1)
plot(npp5[[2]])



### Get nitrate predictors  ----

# DAYS 10 N ----
nit <- rast(paste(re.dir, "Nitrate", "Days_10N", "Days_10N_extended_stack.tif", sep ='/'))
nit

# # crop to NC --
nit2 <- crop(nit, n.extent)
plot(nit2[[1]])

# resample predictors to bathy ----
nit3 <- resample(nit2, depth)

# mask predictors to bathy ----
nit4 <- mask(nit3, depth)
plot(nit4[[1]])



### Get temperature predictors  ----

# Max_Monthly_Anomaly_Upwelling_Temp ----

max_t <- rast(paste(re.dir, "Temperature", "Max_Monthly_Anomaly_Upwelling_Temp", "Max_Monthly_Anomaly_Upwelling_Temp_extended_stack.tif", sep ='/'))
max_t

# # crop to NC --
max_t2 <- crop(max_t, n.extent)
plot(max_t2[[1]])

# resample predictors to bathy ----
max_t3 <- resample(max_t2, depth)

# mask predictors to bathy ----
max_t4 <- mask(max_t3, depth)
plot(max_t4[[1]])



### Get Wave predictors  ----

## UBR_Max ----

ubr.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Environmental/Seagrant_OPC_Environmental_Data/Waves_UBR_30m"

# load raster data --
ubr <- rast(paste(ubr.dir, "UBR_Max_NC-CC_latlon_30mdepth_300res.tif", sep ='/'))
ubr

# # crop to NC --
ubr2 <- crop(ubr, n.extent)
plot(ubr2[[1]])

# resample predictors to bathy ----
ubr3 <- resample(ubr2, depth)

# mask predictors to bathy ----
ubr4 <- mask(ubr3, depth)
plot(ubr4[[1]])


# log ----
ubr5 <- log(ubr4 + 1)
plot(ubr5[[1]])




# make year raster
plot(depth)
year.raster <- classify(depth, cbind(-Inf, Inf, 2004), right=FALSE)
plot(year.raster)
names(year.raster) <- 'survey_year'

year.list <- paste(2004:2021)
length(year.list)



# sites ----

rdf <- as.data.frame(year.raster, xy = T)
head(rdf)
rdf$site <- rdf$y
rdf <- rdf[,-3]
head(rdf)

#site.raster <- rast(rdf, type = 'xyz', crs="EPSG:26910", extent = ext(year2004))

site.raster <- rast(rdf, type = 'xyz', crs="EPSG:4326", extent = ext(year.raster))

site.raster

plot(site.raster)

ext(year.raster)
ext(site.raster)

site.raster2 <- extend(site.raster, year.raster)
plot(site.raster2)

##

# zone ----

sub.dir <- "G:/Shared drives/California Kelp Restoration Project - Seagrant/Data/Susbtrate_aggregated"

plot(depth)
#depth2 <- terra::project(depth, "epsg:4326")


zone.raster <- depth*-1
names(zone.raster) <- 'zone'
plot(zone.raster)

# # # crop to NC --
# zone.raster2 <- crop(zone.raster, n.extent)
# plot(zone.raster2[[1]])
# 
# # resample predictors to bathy ----
# zone.raster3 <- resample(zone.raster2, slope)
# plot(zone.raster3[[1]])
# 
# # mask predictors to bathy ----
# zone.raster4 <- mask(zone.raster3, slope)
# plot(zone.raster4[[1]])

# classify --

rec.m <-  c(-Inf, 10, 1,
            10, Inf, 2)

rclmat <- matrix(rec.m, ncol=3, byrow=TRUE)

zone.raster5 <- classify(zone.raster, rclmat, right=FALSE)
plot(zone.raster5, col = c('red', 'blue'))




##


## LOOP to predict each year using data frame ----

max.urch <- max(dat1$log_den_STRPURAD, na.rm = T) # 8.936298
max.urch

mod.V1 <- bm
summary(mod.V1)

# make list of years --
year.list <- paste(2004:2021)
length(year.list)



# outputs dir ----
# * use an output directory of yours
preds.dir <- paste(cv.dir, "preds_with_predictors", sep ='/')
preds.dir



for (i in 1:length(year.list)) {
  
  # 1. Get env predictors and stack them 
  preds1 <- c(depth, npp.up5[[i+1]], npp5[[i+1]], nit4[[i+6]], max_t4[[i+6]], ubr5[[i]])
  
  # 2. get year and stack it
  year.no <- as.numeric(year.list[i])
  year.r <- classify(year.raster, cbind(-Inf, 0, year.no), right=FALSE)
  
  preds2 <- c(preds1, year.r)
  
  # 3. stack zone
  preds3 <- c(preds2, zone.raster5)
  
  # 4. stack site
  preds4 <- c(preds3, site.raster2)
  
  # 5. name predictors 
  names(preds4) <- c("mean_depth", 
                     "log_Mean_Monthly_NPP_Upwelling", 
                     "log_Min_Monthly_NPP",   
                     "Days_10N",
                     "Max_Monthly_Anomaly_Upwelling_Temp",
                     "log_UBR_Max",
                     "year"     ,
                     "zone",
                     "site_name")
  
  # 6. make data frame for prediction with factors
  df4 <- as.data.frame(preds4, xy = T) %>%
    mutate_at(vars(year, zone, site_name), list(as.factor)) %>%
    mutate(zone = recode_factor(zone, '1' = 'INNER', '2' = 'OUTER')) %>%
    glimpse()
  
  # 7. predict
  year.pred.df <- predict.gam(mod.V1, newdata=df4, type='response', se.fit=T)
  head(year.pred.df)
  
  # preds with predictors
  preds.preds <- df4 %>% 
    data.frame(year.pred.df) %>%
    #dplyr::select(x, y, fit) %>%
    # so urchins do not exceed the maximum recorded
    dplyr::mutate(fit = replace(fit, fit > max.urch, max.urch)) %>%
    na.omit() %>%
    glimpse()
  
  name.csvx <- paste(year.no, "Log_STRPURAD_NC_predictors.csv", sep = '_')
  write.csv(preds.preds, paste(preds.dir, name.csvx, sep = '/'))
  
  # 8. join with df for lats and lons
  preds.all <-  df4 %>% 
    data.frame(year.pred.df) %>%
    dplyr::select(x, y, fit) %>%
    # so urchins do not exceed the maximum recorded
    dplyr::mutate(fit = replace(fit, fit > max.urch, max.urch)) %>%
    glimpse()
  
  # 9. save csv
  name.csv <- paste(year.no, "Log_STRPURAD_NC.csv", sep = '_')
  write.csv(preds.all, paste(preds.dir, name.csv, sep = '/'))
  
  # Save csv without NAs
  preds.all2 <- preds.all %>% na.omit()
  name.csv2 <- paste(year.no, "Log_STRPURAD_NC_noNA.csv", sep = '_')
  write.csv(preds.all2, paste(preds.dir, name.csv2, sep = '/'))
  
  # 10. Rasterize
  crs.p <- "epsg:4326"
  year.prediction <- rast(preds.all, type = 'xyz', crs = crs.p, digits = 6)
  plot(year.prediction)
  
  # 11. save raw raster
  name.raster <- paste(year.no, "Log_STRPURAD_NC.tif", sep = '_')
  writeRaster(year.prediction, paste(preds.dir, name.raster, sep = '/'))
  
}


###

###

### Stack ALL YEAR predictions ----


preds.dir <- paste(o2.dir, "preds.nc", sep ='/')
preds.dir


# list rasters --
n.files <- dir(preds.dir)
# list files in source --
n.files <- list.files(preds.dir, pattern = '.tif', full.names = TRUE)
n.files
length(n.files)
#n.files <- n.files[-c(19, 20)]
# list names to load onto the Environment --
names.list <- list.files(preds.dir, pattern = '.tif')
names.list <- str_replace(names.list, ".tif", "")
length(names.list)
names.list
#names.list <- names.list[-c(19,20)]

# load csv files as a list --
tfiles <- lapply(n.files, rast) # this is a list
tfiles[[1]]
length(tfiles)

# stack them ---
preds.stack <- c()

# use do call to create a raster otherwise it creates a list
preds.stack <- do.call("c", tfiles)
plot(preds.stack[[18]])
plot(preds.stack[[9]])

names(preds.stack) <- paste("log_STRPURAD", paste(2004:2021), sep ='_')
names(preds.stack) 

# # save
# writeRaster(preds.stack, paste(preds.dir, "stack_log_STRPURAD_NC.tif", sep ='/'), overwrite = T)
# 
# ###
# 
# dfs <- as.data.frame(preds.stack, xy = T) %>%
#   glimpse()
# 
# write.csv(dfs, paste(preds.dir, "stack_log_STRPURAD_NC.csv", sep ='/'))

###

plot(preds.stack[[13]])

plot(preds.stack[[18]])

plot(preds.stack[[9]])


